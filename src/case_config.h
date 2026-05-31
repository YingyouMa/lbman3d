#ifndef LBM_AN_CASE_CONFIG_H_
#define LBM_AN_CASE_CONFIG_H_

#include "boundary.h"

#if __has_include("case_config_local.h")
#include "case_config_local.h"
#else
#include "default_case_config.h"
#endif

#endif  // LBM_AN_CASE_CONFIG_H_
