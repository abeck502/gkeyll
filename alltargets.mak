# -*- makefile-gmake -*-
# Generated automatically! Do not edit

.PHONY:
all: moments moments-unit moments-regression ## Build only specified Apps (moments)

.PHONY: install
install: moments-install ## Install gkeyll executable
	cd gkeyll && ${MAKE} -f Makefile-gkeyll install

