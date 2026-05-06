!-----------------------------------------------------------------------------
! (c) Crown copyright 2026 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------
!> @brief   Converts a relative humidity fraction to water vapour mixing ratio.
!> @details At each model level, computes the saturation mixing ratio from the
!!          local potential temperature and Exner pressure, then inverts the
!!          relative humidity expression to obtain the water vapour mixing ratio.
!!
!!          Relative humidity is expected as a fraction where 0 = completely dry
!!          and 1 = saturated.  This is consistent with the definition used in
!!          relative_humidity_kernel_mod.
!!
!!          The inversion formula is derived from the definition used in
!!          relative_humidity_kernel_mod:
!!
!!            rel_hum = (mr_v / mr_sat) * (1 + mr_sat * recip_epsilon)
!!                                       / (1 + mr_v  * recip_epsilon)
!!
!!          Solving for mr_v gives:
!!
!!            mr_v = rel_hum * mr_sat
!!                 / (1 + (1 - rel_hum) * mr_sat * recip_epsilon)

module rh_to_mr_kernel_mod

  use argument_mod,       only: arg_type,         &
                                GH_FIELD,          &
                                GH_SCALAR,         &
                                GH_REAL,           &
                                GH_WRITE, GH_READ, &
                                Wtheta, DOF
  use constants_mod,      only: r_def
  use kernel_mod,         only: kernel_type
  use physics_common_mod, only: qsaturation

  implicit none

  private

  !---------------------------------------------------------------------------
  ! Public types
  !---------------------------------------------------------------------------
  !> The type declaration for the kernel. Contains the metadata needed by the
  !> PSy layer.
  type, public, extends(kernel_type) :: rh_to_mr_kernel_type
    private
    type(arg_type) :: meta_args(7) = (/                               &
         arg_type(GH_FIELD,  GH_REAL, GH_WRITE, Wtheta),             &
         arg_type(GH_FIELD,  GH_REAL, GH_READ,  Wtheta),             &
         arg_type(GH_FIELD,  GH_REAL, GH_READ,  Wtheta),             &
         arg_type(GH_FIELD,  GH_REAL, GH_READ,  Wtheta),             &
         arg_type(GH_SCALAR, GH_REAL, GH_READ),                      &
         arg_type(GH_SCALAR, GH_REAL, GH_READ),                      &
         arg_type(GH_SCALAR, GH_REAL, GH_READ)                       &
         /)
    integer :: operates_on = DOF
  contains
    procedure, nopass :: rh_to_mr_code
  end type rh_to_mr_kernel_type

  public :: rh_to_mr_code

contains

  !> @brief   Converts a relative humidity fraction to water vapour mixing ratio.
  !> @param[in,out] mr_v          Water vapour mixing ratio (kg/kg)
  !> @param[in]     rel_hum       Relative humidity as a fraction (0 = dry, 1 = saturated)
  !> @param[in]     theta         Potential temperature (K)
  !> @param[in]     exner_at_wt   Exner pressure at Wtheta points
  !> @param[in]     p_zero        Reference pressure (Pa)
  !> @param[in]     kappa         Ratio Rd/cp
  !> @param[in]     recip_epsilon Reciprocal of the ratio Rd/Rv (= Rv/Rd)
  subroutine rh_to_mr_code( mr_v, rel_hum, theta, exner_at_wt, &
                             p_zero, kappa, recip_epsilon )

    implicit none

    ! Arguments
    real(kind=r_def), intent(inout) :: mr_v
    real(kind=r_def), intent(in)    :: rel_hum
    real(kind=r_def), intent(in)    :: theta
    real(kind=r_def), intent(in)    :: exner_at_wt
    real(kind=r_def), intent(in)    :: p_zero
    real(kind=r_def), intent(in)    :: kappa
    real(kind=r_def), intent(in)    :: recip_epsilon

    ! Local variables
    real(kind=r_def) :: temperature, pressure, mr_sat

    ! Recover temperature and pressure from potential temperature and Exner.
    temperature = theta * exner_at_wt
    pressure    = p_zero * exner_at_wt ** (1.0_r_def / kappa)

    ! Saturation mixing ratio via Tetens' formula (pressure required in mbar).
    mr_sat = qsaturation(temperature, 0.01_r_def * pressure)

    ! Invert the relative-humidity expression from relative_humidity_kernel_mod:
    !   rel_hum = (mr_v / mr_sat) * (1 + mr_sat * recip_epsilon)
    !                               / (1 + mr_v  * recip_epsilon)
    ! Solved for mr_v:
    !   mr_v = rel_hum * mr_sat / (1 + (1 - rel_hum) * mr_sat * recip_epsilon)
    mr_v = rel_hum * mr_sat / &
           ( 1.0_r_def + (1.0_r_def - rel_hum) * mr_sat * recip_epsilon )

  end subroutine rh_to_mr_code

end module rh_to_mr_kernel_mod
