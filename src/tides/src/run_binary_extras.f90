! ***********************************************************************
!
!   Copyright (C) 2012-2019  Bill Paxton & The MESA Team
!
!   this file is part of mesa.
!
!   mesa is free software; you can redistribute it and/or modify
!   it under the terms of the gnu general library public license as published
!   by the free software foundation; either version 2 of the license, or
!   (at your option) any later version.
!
!   mesa is distributed in the hope that it will be useful, 
!   but without any warranty; without even the implied warranty of
!   merchantability or fitness for a particular purpose.  see the
!   gnu library general public license for more details.
!
!   you should have received a copy of the gnu library general public license
!   along with this software; if not, write to the free software
!   foundation, inc., 59 temple place, suite 330, boston, ma 02111-1307 usa
!
! *********************************************************************** 
      module run_binary_extras 

      use star_lib
      use star_def
      use const_def
      use const_def
      use chem_def
      use num_lib
      use binary_def
      use math_lib
      
      implicit none

      logical, parameter :: dbg = .true.

      ! b% xtra(x_inclination) contains the inclination value (changing at each timestep)
      integer, parameter :: x_inclination = 1
      ! b% xtra(x_idot) contains the derivative of the inclination
      integer, parameter :: x_idot = 2

      ! b% lxtra(lx_convective_envelope) is true when envelope of star is convective
      integer, parameter :: lx_convective_envelope = 1

      ! inclination value
      real(dp) :: initial_inclination
      real(dp), parameter :: min_inclination = 1d-8
      
      contains
      
      subroutine extras_binary_controls(binary_id, ierr)
         integer :: binary_id
         integer, intent(out) :: ierr
         type (binary_info), pointer :: b
         ierr = 0

         call binary_ptr(binary_id, b, ierr)
         if (ierr /= 0) then
            write(*,*) 'failed in binary_ptr'
            return
         end if

         ! apply tides in non-coplanar case
         b% other_sync_spin_to_orbit => sync_non_coplanar
         b% other_edot_tidal => edot_non_coplanar

         ! Set these function pointers to point to the functions you wish to use in
         ! your run_binary_extras. Any which are not set, default to a null_ version
         ! which does nothing.
         b% how_many_extra_binary_history_header_items => how_many_extra_binary_history_header_items
         b% data_for_extra_binary_history_header_items => data_for_extra_binary_history_header_items
         b% how_many_extra_binary_history_columns => how_many_extra_binary_history_columns
         b% data_for_extra_binary_history_columns => data_for_extra_binary_history_columns

         b% extras_binary_startup=> extras_binary_startup
         b% extras_binary_start_step=> extras_binary_start_step
         b% extras_binary_check_model=> extras_binary_check_model
         b% extras_binary_finish_step => extras_binary_finish_step
         b% extras_binary_after_evolve=> extras_binary_after_evolve

         ! Once you have set the function pointers you want, then uncomment this (or set it in your star_job inlist)
         ! to disable the printed warning message,
          b% warn_binary_extra =.false.
         
      end subroutine extras_binary_controls


      subroutine sync_non_coplanar(id, nz, osep, qratio, rl, dt_next, Ftid, sync_type, sync_mode, ierr)
         integer, intent(in) :: id
         integer, intent(in) :: nz
         real(dp), intent(in) :: osep ! orbital separation (cm)
         real(dp), intent(in) :: qratio ! mass_other_star/mass_this_star
         real(dp), intent(in) :: rl ! roche lobe radius (cm)
         real(dp), intent(in) :: dt_next ! next timestep
         real(dp), intent(in) :: Ftid ! efficiency of tidal synchronization. (time scale / Ftid)
         character (len=strlen), intent(in) :: sync_type ! synchronization timescale
         character (len=strlen), intent(in) :: sync_mode ! where to put/take angular momentum
         integer, intent(out) :: ierr
         type(binary_info), pointer :: b
         type(star_info), pointer :: s
         integer :: k
         real(dp) :: moment_of_inertia, rGyr_squared, m, r_phot, porb
         real(dp) :: omega_orb
         real(dp) :: a1, a2, f_sync, t_sync
         real(dp), dimension(nz) :: j_sync, delta_j

         call star_ptr(id, s, ierr)
         if (ierr /= 0) then
            write(*,*) 'failed in star_ptr'
            return
         end if

         call binary_ptr(s% binary_id, b, ierr)
         if (ierr /= 0) then
            write(*,*) 'failed in binary_ptr'
            return
         end if

         s% extra_jdot(:) = 0d0

         if (is_donor(b, s)) then
            m = b% m(b% d_i)
            r_phot = b% r(b% d_i)
         else
            m = b% m(b% a_i)
            r_phot = b% r(b% a_i)
         end if
         porb = b% period
         
         if (is_donor(b, s)) then
            m = b% m(b% d_i)
            r_phot = b% r(b% d_i)
         else
            m = b% m(b% a_i)
            r_phot = b% r(b% a_i)
         end if

         omega_orb = 2d0*pi/b% period
         do k=1,nz
            j_sync(k) = omega_orb*s% i_rot(k)
         end do

         ! compute gyration radius
         moment_of_inertia = dot_product(s% i_rot(:s% nz), s% dm_bar(:s% nz))
         rGyr_squared = (moment_of_inertia/(m*r_phot*r_phot))

         ! compute sync timescale based on different energy transport on the star envelope
         if (b% lxtra(lx_convective_envelope)) then
            t_sync = 3d0 * k_div_T(b, s, .true.) * (qratio*qratio/rGyr_squared) * pow6(r_phot/osep)
         else
            t_sync = 3d0 * k_div_T(b, s, .false.) * (qratio*qratio/rGyr_squared) * pow6(r_phot/osep)
         end if
         t_sync = (1d0 / t_sync) / Ftid

         ! set change in angular momentum on donor star as MESA does
         a1 = f2(b% eccentricity) * cos(b% xtra(x_inclination))
         a2 = pow(1-pow2(b% eccentricity), 1.5d0) * f5(b% eccentricity) * &
            (3d0 + cos(2*b% xtra(x_inclination)))
         do k=1, nz
            delta_j(k) = (1d0 - exp(-a2*dt_next/t_sync))*(s% j_rot(k) - a1/a2*j_sync(k))
         end do

         if (.not. b% doing_first_model_of_run) then
            do k=1,nz
               s% extra_jdot(k) = s% extra_jdot(k) - delta_j(k)/dt_next
            end do
         end if

      end subroutine sync_non_coplanar


      subroutine edot_non_coplanar(binary_id, ierr)
         integer, intent(in) :: binary_id
         integer, intent(out) :: ierr
         type(binary_info), pointer :: b
         type(star_info), pointer :: s
         integer :: k
         real(dp) :: m, porb, r_phot, osep, qratio, omega_s, omega_sync
         real(dp) :: edot

         include 'formats.inc'

         ierr = 0
         call binary_ptr(binary_id, b, ierr)
         if (ierr /= 0) then
            write(*,*) 'failed in binary_ptr'
            return
         end if
         
         s => b% s_donor
         
         porb = b% period
         omega_sync = 2d0*pi/b% period
         omega_s = s% omega_avg_surf
         osep = b% separation
         
         qratio = b% m(b% a_i) / b% m(b% d_i)
         if (is_donor(b, s)) then
            m = b% m(b% d_i)
            r_phot = b% r(b% d_i)
         else
            qratio = 1.0d0/qratio
            m = b% m(b% a_i)
            r_phot = b% r(b% a_i)
         end if

         edot = -27d0 * qratio * (1 + qratio) * pow8(r_phot/osep) * b% eccentricity * &
            pow(1 - pow2(b% eccentricity), -6.5d0) * b% Ftid_1 * &
            k_div_T(b, s, b% lxtra(lx_convective_envelope))
         
         edot = edot * (f3(b% eccentricity) - 11d0/18d0 * omega_s / omega_sync * f4(b% eccentricity) * &
           pow(1-pow2(b% eccentricity),1.5d0) * cos(b% xtra(x_inclination)))

         b% edot_tidal = 0d0 + edot

         if (dbg) write(*,1) 'edot', edot

      end subroutine edot_non_coplanar


      real(dp) function idot_non_coplanar(b, s) result(idot)
         type(binary_info), pointer :: b
         type(star_info), pointer :: s
         real(dp) :: i_step, par
         real(dp) :: porb, omega_sync, omega_s, osep, qratio, m, r_phot
         real(dp) :: rGyr_squared, moment_of_inertia
         
         include 'formats.inc'

         if (b% doing_first_model_of_run .or. b% xtra(x_inclination) <= min_inclination) then
            idot = 0d0
            return
         end if

         i_step = b% xtra(x_inclination)
         
         porb = b% period
         omega_sync = 2d0*pi/b% period
         omega_s = s% omega_avg_surf
         osep = b% separation
         
         qratio = b% m(b% a_i) / b% m(b% d_i)
         if (is_donor(b, s)) then
            m = b% m(b% d_i)
            r_phot = b% r(b% d_i)
         else
            qratio = 1.0d0/qratio
            m = b% m(b% a_i)
            r_phot = b% r(b% a_i)
         end if

         ! calculate the gyration radius squared
         moment_of_inertia = dot_product(s% i_rot(:s% nz), s% dm_bar(:s% nz))
         rGyr_squared = (moment_of_inertia/(m*r_phot*r_phot))

         idot = -3d0 * k_div_T(b, s, b% lxtra(lx_convective_envelope)) * &
            (qratio*qratio / rGyr_squared) * pow6(r_phot/osep) * (omega_sync / omega_s) * &
            pow(1-pow2(b% eccentricity), -6d0) * sin(i_step)

         par = f2(b% eccentricity) - 0.5d0 * f5(b% eccentricity) * &
            ((omega_s/omega_sync) * cos(i_step) * pow(1-pow2(b% eccentricity),1.5d0) + &
            r_phot*r_phot * osep * omega_s*omega_s * rGyr_squared * (1 - pow2(b% eccentricity)) / &
            (b% m(b% a_i) * standard_cgrav))

         idot = par * idot

         if (dbg) write(*,1) 'idot', idot

      end function idot_non_coplanar


      logical function is_convective(s)
         type(star_info), pointer :: s
         real(dp) :: m_conv, f_conv
         integer :: k

         include 'formats.inc'

         is_convective = .false.
         
         if (s% he_core_mass > 0d0 .and. s% center_h1 < 1d-6) then
            if (s% he_core_k == 1 .or. (s% star_mass - s% he_core_mass) < 1d-3) then
               m_conv = 0d0
               f_conv = 0d0
            else
               do k=1,s% he_core_k
                  if (s% mixing_type(k) == convective_mixing) m_conv = m_conv + s% dm(k)
               end do
               f_conv = (m_conv / Msun) / (s% star_mass - s% he_core_mass)
            end if
         else
            m_conv = 0d0
            f_conv = 0d0
         end if
         
         if (f_conv > 0.2) is_convective = .true.

         if (dbg) write(*,14) 'is_convective', is_convective

      end function is_convective


      real(dp) function k_div_T(b, s, has_convective_envelope)
         type(binary_info), pointer :: b
         type(star_info), pointer :: s
         logical, intent(in) :: has_convective_envelope
         integer :: k
         real(dp) osep, qratio, m, r_phot,porb, m_env, r_env, tau_conv, P_tid, f_conv

         ! k/T computed as in Hurley, J., Tout, C., Pols, O. 2002, MNRAS, 329, 897
         ! Kudos to Francesca Valsecchi for help implementing and testing this

         k_div_T = 0d0

         osep = b% separation
         qratio = b% m(b% a_i) / b% m(b% d_i)
         if (is_donor(b, s)) then
            m = b% m(b% d_i)
            r_phot = b% r(b% d_i)
         else
            qratio = 1d0/qratio
            m = b% m(b% a_i)
            r_phot = b% r(b% a_i)
         end if
         porb = b% period

         if (has_convective_envelope) then
            m_env = 0d0
            r_env = 0d0
            do k=1, s% nz
               if (s% mixing_type(k) /= convective_mixing .and. &
                   s% rho(k) > 1d5*s% rho(1)) then
                  r_env = (r_phot - s% r(k))/Rsun
                  m_env = (s% m(1) - s% m(k))/Msun
                  exit
               end if
            end do
            tau_conv = 0.431d0*pow(m_env*r_env* &
               (r_phot/Rsun-r_env/2d0)/3d0/s% L_phot,one_third) * secyer
            if (1d0/porb /= s% omega_avg_surf/(2d0*pi)) then
               P_tid = 1d0/abs(1d0/porb-s% omega_avg_surf/(2d0*pi))
               f_conv = min(1.0d0, pow(P_tid/(2d0*tau_conv),b% tidal_reduction))
            else
               f_conv = 1d0
            end if

            k_div_T = 2d0/21d0*f_conv/tau_conv*m_env/(m/Msun)
         else
            !NOTE:There is a typo in eq. (42) of Hurley+ 2002,
            !correct expression is given in footnote 3 of
            !Sepinsky+ 2007
            k_div_T = 1.9782d4*sqrt(m*r_phot*r_phot/pow5(osep)/(Msun/pow3(Rsun)))
            k_div_T = k_div_T*pow(1d0+qratio,5d0/6d0)
            k_div_T = k_div_T*1.592d-9*pow(m/Msun,2.84d0)/secyer
         end if

      end function k_div_T


      real(dp) function f2(e)
         real(dp), intent(in) :: e

         f2 = 1d0

         ! Hut 1981, A&A, 99, 126, definition of f2 after eq. 11
         if (e > 0d0) then
             f2 = 1d0 + 15d0/2d0*pow2(e) + 45d0/8d0*pow4(e) + 5d0/16d0*pow6(e)
         end if

      end function f2


      real(dp) function f3(e)
         real(dp), intent(in) :: e

         f3 = 1d0

         ! Hut 1981, A&A, 99, 126, definition of f3 after eq. 11
         if (e > 0d0) then
             f3 = 1d0 + 15d0/4d0*pow2(e) + 15d0/8d0*pow4(e) + 5d0/64d0*pow6(e)
         end if

      end function f3


      real(dp) function f4(e)
         real(dp), intent(in) :: e

         f4 = 1d0

         ! Hut 1981, A&A, 99, 126, definition of f4 after eq. 11
         if (e > 0d0) then
             f4 = 1d0 + 3d0/2d0*pow2(e) + 1d0/8d0*pow4(e)
         end if

      end function f4


      real(dp) function f5(e)
         real(dp), intent(in) :: e

         f5 = 1d0

         ! Hut 1981, A&A, 99, 126, definition of f5 after eq. 11
         if (e > 0d0) then
             f5 = 1d0 + 3d0*pow2(e) + 3d0/8d0*pow4(e)
         end if

      end function f5


      integer function how_many_extra_binary_history_header_items(binary_id)
         use binary_def, only: binary_info
         integer, intent(in) :: binary_id
         how_many_extra_binary_history_header_items = 0
      end function how_many_extra_binary_history_header_items


      subroutine data_for_extra_binary_history_header_items( &
           binary_id, n, names, vals, ierr)
         type (binary_info), pointer :: b
         integer, intent(in) :: binary_id, n
         character (len=maxlen_binary_history_column_name) :: names(n)
         real(dp) :: vals(n)
         integer, intent(out) :: ierr
         ierr = 0
         call binary_ptr(binary_id, b, ierr)
         if (ierr /= 0) then
            write(*,*) 'failed in binary_ptr'
            return
         end if
      end subroutine data_for_extra_binary_history_header_items


      integer function how_many_extra_binary_history_columns(binary_id)
         use binary_def, only: binary_info
         integer, intent(in) :: binary_id

         how_many_extra_binary_history_columns = 2

      end function how_many_extra_binary_history_columns


      subroutine data_for_extra_binary_history_columns(binary_id, n, names, vals, ierr)
         type (binary_info), pointer :: b
         integer, intent(in) :: binary_id
         integer, intent(in) :: n
         character (len=maxlen_binary_history_column_name) :: names(n)
         real(dp) :: vals(n)
         integer, intent(out) :: ierr

         ierr = 0
         call binary_ptr(binary_id, b, ierr)
         if (ierr /= 0) then
            write(*,*) 'failed in binary_ptr'
            return
         end if

         names(1) = 'inclination'
         names(2) = 'idot'

         vals(1) = b% xtra(x_inclination) * rad2a
         vals(2) = b% xtra(x_idot) * rad2a

      end subroutine data_for_extra_binary_history_columns
      
      
      integer function extras_binary_startup(binary_id,restart,ierr)
         type (binary_info), pointer :: b
         integer, intent(in) :: binary_id
         integer, intent(out) :: ierr
         logical, intent(in) :: restart

         call binary_ptr(binary_id, b, ierr)
         if (ierr /= 0) then ! failure in  binary_ptr
            return
         end if
         
         ! for now, just use it as input
         initial_inclination = b% s1% x_ctrl(1)

         if (.not. restart) then
            b% lxtra(lx_convective_envelope) = .false.
            b% xtra(x_inclination) = initial_inclination * a2rad
            b% xtra(x_idot) = 0d0
         end if

         extras_binary_startup = keep_going

      end function  extras_binary_startup
     

      integer function extras_binary_start_step(binary_id,ierr)
         type (binary_info), pointer :: b
         integer, intent(in) :: binary_id
         integer, intent(out) :: ierr
         type(star_info), pointer :: s
         real(dp) :: idot

         include 'formats.inc'

         call binary_ptr(binary_id, b, ierr)
         if (ierr /= 0) then ! failure in  binary_ptr
            return
         end if
         
         extras_binary_start_step = keep_going
         
         ! check if envelope is convective or not
         s => b% s_donor
         b% lxtra(lx_convective_envelope) = is_convective(s)

         ! to update inclination based on Repetto & Nelemans (2014)
         ! first we compute di
         idot = idot_non_coplanar(b, s)

         ! then we can update inclination
         b% xtra(x_inclination) = b% xtra(x_inclination) + idot * b% time_step * secyer
         b% xtra(x_idot) = idot

         ! do not go beyond minimum inclination value
         if (b% xtra(x_inclination) < min_inclination) b% xtra(x_inclination) = min_inclination

         if (dbg) write(*,1) 'inclination', b% xtra(x_inclination) * rad2a

      end function extras_binary_start_step
      

      !Return either keep_going, retry or terminate
      integer function extras_binary_check_model(binary_id)
         type (binary_info), pointer :: b
         integer, intent(in) :: binary_id
         integer :: ierr

         call binary_ptr(binary_id, b, ierr)
         if (ierr /= 0) then ! failure in  binary_ptr
            return
         end if

         extras_binary_check_model = keep_going
        
      end function extras_binary_check_model
      
      
      ! returns either keep_going or terminate.
      ! note: cannot request retry; extras_check_model can do that.
      integer function extras_binary_finish_step(binary_id)
         type (binary_info), pointer :: b
         integer, intent(in) :: binary_id
         integer :: ierr

         call binary_ptr(binary_id, b, ierr)
         if (ierr /= 0) then ! failure in  binary_ptr
            return
         end if

         extras_binary_finish_step = keep_going
         
      end function extras_binary_finish_step
     

      subroutine extras_binary_after_evolve(binary_id, ierr)
         type (binary_info), pointer :: b
         integer, intent(in) :: binary_id
         integer, intent(out) :: ierr

         call binary_ptr(binary_id, b, ierr)
         if (ierr /= 0) then ! failure in  binary_ptr
            return
         end if 
         
      end subroutine extras_binary_after_evolve     
      

      end module run_binary_extras
