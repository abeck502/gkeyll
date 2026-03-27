#pragma once

#include <math.h>

#include <gkyl_array.h>
#include <gkyl_eqn_type.h>
#include <gkyl_range.h>
#include <gkyl_rect_grid.h>
#include <gkyl_util.h>

// Identifiers for different viscosity types.
enum gkyl_viscosity_type
{
  GKYL_VISC_NONE = 0,
  GKYL_VISC_CONST = 1 << 0
};

struct gkyl_moment_viscosity_data {
  enum gkyl_eqn_type type_eqn; // equation type
  enum gkyl_viscosity_type type_visc; // which viscosity terms are enabled
  double mu; // dynamic viscosity coefficient
};

struct gkyl_moment_viscosity_inp {
  const struct gkyl_rect_grid *grid; // grid on which to solve equations
  int nfluids; // number of fluids
  struct gkyl_moment_viscosity_data param[GKYL_MAX_SPECIES]; // species data
};

// Object type
typedef struct gkyl_moment_viscosity gkyl_moment_viscosity;

/**
 * Create new updater to update fluid equations with viscosity
 * transport terms.
 * Returns RHS for accumulation in a forward Euler method.
 *
 * @param inp Input parameters to updater
 */
gkyl_moment_viscosity* gkyl_moment_viscosity_new(struct gkyl_moment_viscosity_inp inp);

/**
 * Compute RHS contribution from viscosity transport terms in the
 * multi-fluid system. The update_rng MUST be a sub-range of the
 * range on which the array is defined. That is, it must be either
 * the same range as the array range, or one created using the
 * gkyl_sub_range_init method.
 *
 * @param ves New viscosity updater object.
 * @param visc_vars_rng Range on which to compute viscosity variables (cell nodes)
 * @param update_rng Range on which to compute update.
 * @param fluid Input array of fluid variables (array size: nfluids)
 * @param cflrate CFL scalar rate (frequency: units of 1/[T]) (array size: nfluids)
 * @param visc_vars Array for storing intermediate computation of viscosity variables (cell nodes)
 * @param rhs RHS output (NOTE: Returns RHS output of all nfluids)
 */
void gkyl_moment_viscosity_advance(const gkyl_moment_viscosity *ves,
  struct gkyl_range visc_vars_range, struct gkyl_range update_range,
  struct gkyl_array *fluid[GKYL_MAX_SPECIES],
  struct gkyl_array *cflrate[GKYL_MAX_SPECIES],
  struct gkyl_array *visc_vars[GKYL_MAX_SPECIES], struct gkyl_array *rhs[GKYL_MAX_SPECIES]);

/**
 * Delete updater.
 *
 * @param ves Updater to delete.
 */
void gkyl_moment_viscosity_release(gkyl_moment_viscosity *ves);
