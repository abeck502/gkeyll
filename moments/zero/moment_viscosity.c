#include <gkyl_alloc.h>
#include <gkyl_array_ops.h>
#include <gkyl_moment_viscosity.h>
#include <gkyl_moment_non_ideal_priv.h>
#include <gkyl_wv_euler_priv.h>

// 1D stencil locations (L: lower, U: upper)
enum loc_1d {
  L_1D, U_1D
};

// 2D stencil locations (L: lower, U: upper)
enum loc_2d {
  LL_2D, LU_2D,
  UL_2D, UU_2D
};

// 3D stencil locations (L: lower, U: upper)
enum loc_3d {
  LLL_3D, LLU_3D,
  LUL_3D, LUU_3D,
  ULL_3D, ULU_3D,
  UUL_3D, UUU_3D
};

struct gkyl_moment_viscosity {
  struct gkyl_rect_grid grid; // grid object
  int ndim; // number of dimensions
  int nfluids; // number of fluids in multi-fluid system
  struct gkyl_moment_viscosity_data param[GKYL_MAX_SPECIES]; // struct of fluid parameters
};

static void
create_offsets_vertices(const struct gkyl_range *range, long offsets[])
{
  // box spanning stencil
  struct gkyl_range box3;
  gkyl_range_init(&box3, range->ndim, (int[]) { -1, -1, -1 }, (int[]) { 0, 0, 0 });

  struct gkyl_range_iter iter3;
  gkyl_range_iter_init(&iter3, &box3);

  // construct list of offsets
  int count = 0;
  while (gkyl_range_iter_next(&iter3))
    offsets[count++] = gkyl_range_offset(range, iter3.idx);
}

static void
create_offsets_centers(const struct gkyl_range *range, long offsets[])
{
  // box spanning stencil
  struct gkyl_range box3;
  gkyl_range_init(&box3, range->ndim, (int[]) { 0, 0, 0 }, (int[]) { 1, 1, 1 });

  struct gkyl_range_iter iter3;
  gkyl_range_iter_init(&iter3, &box3);

  // construct list of offsets
  int count = 0;
  while (gkyl_range_iter_next(&iter3))
    offsets[count++] = gkyl_range_offset(range, iter3.idx);
}

static bool
has_viscosity(const gkyl_moment_viscosity *ves, int n)
{
  return (ves->param[n].type_visc & GKYL_VISC_CONST) && (ves->param[n].mu > 0.0);
}

static void
visc_calc_vars(const gkyl_moment_viscosity *ves,
  const double *fluid_d[][GKYL_MAX_SPECIES],
  double *cflrate[GKYL_MAX_SPECIES], double *visc_d[GKYL_MAX_SPECIES])
{
  const int nfluids = ves->nfluids;
  const int ndim = ves->ndim;

  (void) cflrate;

  for (int n = 0; n < nfluids; ++n)
    for (int k = 0; k <= (int) QZ; ++k)
      visc_d[n][k] = 0.0;

  double eta_avg[GKYL_MAX_SPECIES] = {};
  double w[GKYL_MAX_SPECIES][6] = {};
  double u_avg[GKYL_MAX_SPECIES][3] = {};

  if (ndim == 1) {
    double dx = ves->grid.dx[0];

    // Derived quantities
    double u[2][GKYL_MAX_SPECIES][3] = {}; // Flow for each species, ux, uy, & uz

    for (int j = L_1D; j <= U_1D; ++j) {
      for (int n = 0; n < nfluids; ++n) {
        double rho = fluid_d[j][n][RHO];
        if (rho != 0.0) {
          for (int k = 0; k < 3; ++k)
            u[j][n][k] = fluid_d[j][n][MX + k] / rho;
        }
        else {
          u[j][n][0] = 0.0;
          u[j][n][1] = 0.0;
          u[j][n][2] = 0.0;
        }
      }
    }

    for (int n = 0; n < nfluids; ++n) {
      if (!has_viscosity(ves, n))
        continue;

      // Viscosity coefficient at cell edges (using harmonic average)
      eta_avg[n] = calc_harmonic_avg_1D(ves->param[n].mu, ves->param[n].mu);

      // Rate of strain tensor at cell edges
      calc_ros_1D(dx, u[L_1D][n], u[U_1D][n], w[n]);

      // Calculate viscous heating if energy variable exists
      if (ves->param[n].type_eqn == GKYL_EQN_EULER)
        for (int k = 0; k < 3; ++k)
          u_avg[n][k] = calc_arithm_avg_1D(u[L_1D][n][k], u[U_1D][n][k]);
    }
  }
  else if (ndim == 2) {
    double dx = ves->grid.dx[0];
    double dy = ves->grid.dx[1];

    // Derived quantities
    double u[4][GKYL_MAX_SPECIES][3] = {}; // Flow for each species, ux, uy, & uz

    for (int j = LL_2D; j <= UU_2D; ++j) {
      for (int n = 0; n < nfluids; ++n) {
        double rho = fluid_d[j][n][RHO];
        if (rho != 0.0) {
          for (int k = 0; k < 3; ++k)
            u[j][n][k] = fluid_d[j][n][MX + k] / rho;
        }
        else {
          u[j][n][0] = 0.0;
          u[j][n][1] = 0.0;
          u[j][n][2] = 0.0;
        }
      }
    }

    for (int n = 0; n < nfluids; ++n) {
      if (!has_viscosity(ves, n))
        continue;

      // Viscosity coefficient at cell vertices (using harmonic average)
      eta_avg[n] = calc_harmonic_avg_2D(ves->param[n].mu, ves->param[n].mu, ves->param[n].mu, ves->param[n].mu);

      // Rate of strain tensor at cell vertices
      calc_ros_2D(dx, dy, u[LL_2D][n], u[LU_2D][n], u[UL_2D][n], u[UU_2D][n], w[n]);

      // Calculate viscous heating if energy variable exists
      if (ves->param[n].type_eqn == GKYL_EQN_EULER)
        for (int k = 0; k < 3; ++k)
          u_avg[n][k] = calc_arithm_avg_2D(u[LL_2D][n][k], u[LU_2D][n][k], u[UL_2D][n][k], u[UU_2D][n][k]);
    }
  }

  for (int n = 0; n < nfluids; ++n) {
    if (!has_viscosity(ves, n))
      continue;

    // Total viscous stress tensor
    for (int k = 0; k < 6; ++k)
      visc_d[n][PIXX + k] = -eta_avg[n] * w[n][k];

    // Pi dot u
    double Piu[3] = { visc_d[n][PIXX]*u_avg[n][0] + visc_d[n][PIXY]*u_avg[n][1] + visc_d[n][PIXZ]*u_avg[n][2],
                      visc_d[n][PIXY]*u_avg[n][0] + visc_d[n][PIYY]*u_avg[n][1] + visc_d[n][PIYZ]*u_avg[n][2],
                      visc_d[n][PIXZ]*u_avg[n][0] + visc_d[n][PIYZ]*u_avg[n][1] + visc_d[n][PIZZ]*u_avg[n][2] };

    // Total viscous heating
    if (ves->param[n].type_eqn == GKYL_EQN_EULER)
      for (int k = 0; k < 3; ++k)
        visc_d[n][QX + k] = Piu[k];
  }
}

static void
visc_calc_update(const gkyl_moment_viscosity *ves,
  const double *visc_d[][GKYL_MAX_SPECIES], double *rhs[GKYL_MAX_SPECIES])
{
  int nfluids = ves->nfluids;
  const int ndim = ves->ndim;
  double div_pi[GKYL_MAX_SPECIES][3] = {};
  double div_q[GKYL_MAX_SPECIES] = {};

  if (ndim == 1) {
    const double dx = ves->grid.dx[0];
    double pi[2][GKYL_MAX_SPECIES][6] = {};
    double q[2][GKYL_MAX_SPECIES][3] = {};

    for (int n = 0; n < nfluids; ++n) {
      for (int j = L_1D; j <= U_1D; ++j)
        for (int k = 0; k < 6; ++k)
          pi[j][n][k] = visc_d[j][n][PIXX + k];

      div_pi[n][0] = calc_sym_grad_1D(dx, pi[L_1D][n][0], pi[U_1D][n][0]);
      div_pi[n][1] = calc_sym_grad_1D(dx, pi[L_1D][n][1], pi[U_1D][n][1]);
      div_pi[n][2] = calc_sym_grad_1D(dx, pi[L_1D][n][2], pi[U_1D][n][2]);

      rhs[n][RHO] = 0.0;
      rhs[n][MX] = -div_pi[n][0];
      rhs[n][MY] = -div_pi[n][1];
      rhs[n][MZ] = -div_pi[n][2];

      // If energy variable exists, increment viscous heating
      if (ves->param[n].type_eqn == GKYL_EQN_EULER) {
        for (int j = L_1D; j <= U_1D; ++j)
          for (int k = 0; k < 3; ++k)
            q[j][n][k] = visc_d[j][n][QX + k];

        div_q[n] = calc_sym_grad_1D(dx, q[L_1D][n][0], q[U_1D][n][0]);
        rhs[n][ER] = -div_q[n];
      }
    }
  }
  else if (ndim == 2) {
    const double dx = ves->grid.dx[0];
    const double dy = ves->grid.dx[1];
    double pi[4][GKYL_MAX_SPECIES][6] = {};
    double q[4][GKYL_MAX_SPECIES][3] = {};

    for (int n = 0; n < nfluids; ++n) {
      for (int j = LL_2D; j <= UU_2D; ++j)
        for (int k = 0; k < 6; ++k)
          pi[j][n][k] = visc_d[j][n][PIXX + k];

      div_pi[n][0] = calc_sym_gradx_2D(dx, pi[LL_2D][n][0], pi[LU_2D][n][0], pi[UL_2D][n][0], pi[UU_2D][n][0])
                     + calc_sym_grady_2D(dy, pi[LL_2D][n][1], pi[LU_2D][n][1], pi[UL_2D][n][1], pi[UU_2D][n][1]);

      div_pi[n][1] = calc_sym_gradx_2D(dx, pi[LL_2D][n][1], pi[LU_2D][n][1], pi[UL_2D][n][1], pi[UU_2D][n][1])
                     + calc_sym_grady_2D(dy, pi[LL_2D][n][3], pi[LU_2D][n][3], pi[UL_2D][n][3], pi[UU_2D][n][3]);

      div_pi[n][2] = calc_sym_gradx_2D(dx, pi[LL_2D][n][2], pi[LU_2D][n][2], pi[UL_2D][n][2], pi[UU_2D][n][2])
                     + calc_sym_grady_2D(dy, pi[LL_2D][n][4], pi[LU_2D][n][4], pi[UL_2D][n][4], pi[UU_2D][n][4]);

      rhs[n][RHO] = 0.0;
      rhs[n][MX] = -div_pi[n][0];
      rhs[n][MY] = -div_pi[n][1];
      rhs[n][MZ] = -div_pi[n][2];

      // If energy variable exists, increment viscous heating
      if (ves->param[n].type_eqn == GKYL_EQN_EULER) {
        for (int j = LL_2D; j <= UU_2D; ++j)
          for (int k = 0; k < 3; ++k)
            q[j][n][k] = visc_d[j][n][QX + k];

        div_q[n] = calc_sym_gradx_2D(dx, q[LL_2D][n][0], q[LU_2D][n][0], q[UL_2D][n][0], q[UU_2D][n][0])
                   + calc_sym_grady_2D(dy, q[LL_2D][n][1], q[LU_2D][n][1], q[UL_2D][n][1], q[UU_2D][n][1]);
        rhs[n][ER] = -div_q[n];
      }
    }
  }
}

gkyl_moment_viscosity*
gkyl_moment_viscosity_new(struct gkyl_moment_viscosity_inp inp)
{
  gkyl_moment_viscosity *up = gkyl_malloc(sizeof(gkyl_moment_viscosity));

  up->grid = *(inp.grid);
  up->ndim = up->grid.ndim;
  up->nfluids = inp.nfluids;
  for (int n = 0; n < inp.nfluids; ++n)
    up->param[n] = inp.param[n];

  return up;
}

void
gkyl_moment_viscosity_advance(const gkyl_moment_viscosity *ves,
  struct gkyl_range visc_vars_range, struct gkyl_range update_range,
  struct gkyl_array *fluid[GKYL_MAX_SPECIES],
  struct gkyl_array *cflrate[GKYL_MAX_SPECIES],
  struct gkyl_array *visc_vars[GKYL_MAX_SPECIES], struct gkyl_array *rhs[GKYL_MAX_SPECIES])
{
  int nfluids = ves->nfluids;
  int ndim = update_range.ndim;
  long sz[] = { 2, 4, 8 };

  long offsets_vertices[sz[ndim-1]];
  create_offsets_vertices(&visc_vars_range, offsets_vertices);

  long offsets_centers[sz[ndim-1]];
  create_offsets_centers(&update_range, offsets_centers);

  const double* fluid_d[sz[ndim-1]][GKYL_MAX_SPECIES];
  double *cflrate_d[GKYL_MAX_SPECIES];
  double *visc_vars_d[GKYL_MAX_SPECIES];
  const double* visc_vars_up[sz[ndim-1]][GKYL_MAX_SPECIES];
  double *rhs_d[GKYL_MAX_SPECIES];

  struct gkyl_range_iter iter_vertex;
  gkyl_range_iter_init(&iter_vertex, &visc_vars_range);
  while (gkyl_range_iter_next(&iter_vertex)) {

    long linc_vertex = gkyl_range_idx(&visc_vars_range, iter_vertex.idx);
    long linc_center = gkyl_range_idx(&update_range, iter_vertex.idx);

    for (int i = 0; i < sz[ndim-1]; ++i)
      for (int n = 0; n < nfluids; ++n)
        fluid_d[i][n] = gkyl_array_cfetch(fluid[n], linc_center + offsets_vertices[i]);

    for (int n = 0; n < nfluids; ++n) {
      visc_vars_d[n] = gkyl_array_fetch(visc_vars[n], linc_vertex);
      cflrate_d[n] = gkyl_array_fetch(cflrate[n], linc_center);
    }

    visc_calc_vars(ves, fluid_d, cflrate_d, visc_vars_d);
  }

  struct gkyl_range_iter iter_center;
  gkyl_range_iter_init(&iter_center, &update_range);
  while (gkyl_range_iter_next(&iter_center)) {

    long linc_vertex = gkyl_range_idx(&visc_vars_range, iter_center.idx);
    long linc_center = gkyl_range_idx(&update_range, iter_center.idx);

    for (int i = 0; i < sz[ndim-1]; ++i)
      for (int n = 0; n < nfluids; ++n)
        visc_vars_up[i][n] = gkyl_array_fetch(visc_vars[n], linc_vertex + offsets_centers[i]);

    for (int n = 0; n < nfluids; ++n)
      rhs_d[n] = gkyl_array_fetch(rhs[n], linc_center);

    visc_calc_update(ves, visc_vars_up, rhs_d);
  }
}

void
gkyl_moment_viscosity_release(gkyl_moment_viscosity* up)
{
  free(up);
}
