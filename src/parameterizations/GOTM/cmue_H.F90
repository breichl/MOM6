#include"cppdefs.h"
!-----------------------------------------------------------------------
!BOP
!
! !ROUTINE: The non-local, exact weak-equilibrium stability function \label{sec:cmueA}
!
! !INTERFACE:
   subroutine cmue_h(nlev)

! !DESCRIPTION:
!
!  The solution of \eq{bijVertical} and \eq{giVertical} has the shape indicated
!  by \eq{b13}. This subroutine is used to update the quantities
!  $c_\mu$, $c'_\mu$ and $\Gamma$, defined in \eq{b13}, from which all turbulent
!  fluxes can be computed. The non-linear terms ${\cal N}$ and ${\cal N}_b$ are updated
!  by evaluating the right hand side of \eq{NandNb} at the old time step.
!
!  The numerators and the denominator appearing in \eq{cm}
!  are polynomials of the form
!  \begin{equation}
!   \label{vdng}
!    \begin{array}{rcl}
!    D         &=& d_0
!               +  d_1 \overline{N}^2  + d_2 \overline{S}^2
!               +  d_3 \overline{N}^2 \overline{S}^2
!               + d_4 \overline{N}^4   + d_5 \overline{S}^4      \comma \\[3mm]
!    N_n        &=& n_0
!               +  n_1 \overline{N}^2  + n_2 \overline{S}^2
!               +  n_3 \overline{T}                              \comma \\[3mm]
!    N_b       &=& n_{b0}
!               +  n_{b1} \overline{N}^2 + n_{b2} \overline{S}^2 \comma \\[3mm]
!    N_\Gamma  &=& ( g_0
!               +  g_1 \overline{N}^2  + g_2 \overline{S}^2 ) \overline{T}
!   \point
!   \end{array}
!  \end{equation}
!
!  The coefficients of $D$ are given by
!  \begin{equation}
!   \label{vdi}
!   \begin{array}{rcl}
!     d_0 &=& 36 {\cal N}^3 {\cal N}_b^2                           \comma \\[3mm]
!     d_1 &=& 84 a_5 a_{b3} {\cal N}^2 {\cal N}_b                  \comma \\[3mm]
!     d_2 &=&  9 (a_{b2}^2 - a_{b1}^2) {\cal N}^3
!           + 12 ( 3 a_3^2 - a_2^2) {\cal N}   {\cal N}_b^2        \comma \\[3mm]
!     d_3 &=& 12 (a_2 a_{b1} - 3 a_3 a_{b2} ) a_5 a_{b3} {\cal N}
!           + 12 ( a_3^2 - a_2^2)  a_5 a_{b3} {\cal N}_b           \comma \\[3mm]
!     d_4 &=& 48 a_5^2 a_{b3}^2 {\cal N}                           \comma \\[3mm]
!     d_5 &=&  3 ( 3 a_3^2 - a_2^2)(a_{b2}^2 - a_{b1}^2) {\cal N}
!     \point
!   \end{array}
!  \end{equation}
!  The coefficients of the numerators $N_n$ and $N_b$ can be expressed as
!  \begin{equation}
!   \label{vni}
!   \begin{array}{rcl}
!     n_0 &=& 36 a_1 {\cal N}^2 {\cal N}_b^2                       \comma \\[3mm]
!     n_1 &=&-12 a_5 a_{b3} (a_{b1}+a_{b2}) {\cal N}^2
!           -  8 a_5 a_{b3} (-6 a_1+a_2+3 a_3) {\cal N} {\cal N}_b \comma \\[3mm]
!     n_2 &=& 9 a_1 (a_{b2}^2 - a_{b1}^2){\cal N}^2                \comma \\[3mm]
!     n_3 &=& 36 a_5 a_{b4} (a_{b1}+a_{b2}) {\cal N}^2
!           + 24 a_5 a_{b4} (a_2+3 a_3) {\cal N} {\cal N}_b        \comma
!   \end{array}
!  \end{equation}
!  \begin{equation}
!   \label{vnbi}
!   \begin{array}{rcl}
!     n_{b0} &=& 12 a_{b3} {\cal N}^3 {\cal N}_b                   \comma \\[3mm]
!     n_{b1} &=& 12 a_5 a_{b3}^2 {\cal N}^2                        \comma \\[3mm]
!     n_{b2} &=&  9 a_1 a_{b3} (a_{b1}-a_{b2}) {\cal N}^2
!              +    a_{b3} (6 a_1 (a_2-3 a_3)
!              -    4 (a_2^2-3 a_3^2) ) {\cal N} {\cal N}_b        \comma
!   \end{array}
!  \end{equation}
!  and the numerator of the term $\Gamma$ is
!  \begin{equation}
!   \label{vgi}
!   \begin{array}{rcl}
!     g_0 &=& 36 a_{b4} {\cal N}^3 {\cal N}_b                      \comma \\[3mm]
!     g_1 &=& 36 a_5 a_{b3}  a_{b4} {\cal N}^2                     \comma \\[3mm]
!     g_2 &=& 12 a_{b4} ( 3 a_3^2 - a_2^2) {\cal N}  {\cal N}_b    \point
!   \end{array}
!  \end{equation}

!
! !USES:
   use turbulence, only: an,as,at,av,astk
   use turbulence, only: cmue1,cmue2,cmue3

   IMPLICIT NONE
!
! !INPUT PARAMETERS:
!
!  number of vertical layers
   integer, intent(in)       :: nlev
!
! !BUGS:
! Test stage. Do not yet use.
!
! !REVISION HISTORY:
!  Original author(s): Lars Umlauf
!
!EOP
!-----------------------------------------------------------------------!
! !LOCAL VARIABLES:
     integer                 ::   i
     REALTYPE                ::   A1, A2
     REALTYPE                ::   B1, B2
     REALTYPE                ::   C1, C2, C3
     REALTYPE                ::   CS1, CS2
     REALTYPE                ::   DS0, DS1, DS2
     REALTYPE                ::   DM0, DM1, DM2, DM3
     REALTYPE                ::   DH0, DH1, DH2, DH3
     REALTYPE                ::   SHup, SHdn, SH, SM, SMS
!
!-----------------------------------------------------------------------
!BOC

     A1 = 0.92
     A2 = 0.74
     B1 = 16.6
     B2 = 10.1
     C1 = 0.08
     C2 = 0.7
     C3 = 0.2
     CS1 = 0.0!C1
     CS2 = 0.0!C2
     
     DS0 = A1*(1.-6.*A1/B1-3.*CS1)
     DM0 = A1*(1.-6.*A1/B1-3.*C1)
     DH0 = A2*(1.-6.*A1/B1)
     cmue1(:)=0.0;cmue2(:)=0.0;cmue3(:)=0.0
     DO I = 1, NLEV-1

!        print*,'-',AS(I),AN(I),AV(I)
        AN(I)=min(AN(I),0.015)
        AN(I)=max(-0.28,AN(I))
        AV(I)=min(AV(I),0.015)
        ASTK(I)=min(ASTK(I),0.015)

        DS1 = 1. - 9.*A1 * (A2*AN(I) + A1*AV(I))
        DS2 = -9.*A1*A2*CS2*AN(I)
        DM1 = 1. - 9.*A1*(A2*AN(I)+4.*A1*AV(I))
        DM2 = 9.*A1*(2.*A1+A2*(1.-C2))*AN(I)
        DM3 = 27.*A1*A1*ASTK(I)
        DH1 = 1.-3.*A2*((6.*A1+B2*(1.-C3))*AN(I) + 3.*A2*(1.-C2)*AV(I) - 3.*A2*CS2*ASTK(I))
        DH2 = 9.*A2*AV(I)*(2.*A1+A2)
        DH3 = 9.*A2*ASTK(I)*(2.*A1+A2)

        SHup = DH0*DM1*DS1 + DH2*DM0*DS1 + (DH3*DM1 + DM3*DH2)*DS0
        SHdn = DH1*DM1*DS1 - DH2*DM2*DS1 - (DH3*DM1 + DM3*DH2)*DS2
        SH = SHup/SHdn
        SMS = (DS0 + DS2*SH) / DS1
        SM = (DM0 + DM2*SH + DM3*SMS) / DM1
        cmue1(I) = max(SM,1.e-10)
        cmue2(I) = max(SH,1.e-10)
        cmue3(I) = max(SMS,1.e-10)
   
        cmue3(I) = a2*((9.*a1*a2*an(i)-1.)*(1.-6*a1/b1))  / &
                      ((9.*a1*a2*an(i)-1.)*(1.-3*a2*an(i)*(6*a1+b2*(1-c3))))

     end do

     

     !print*,'h--------------'
     !print*,cmue1(NLEV-5:NLEV)
     !print*,cmue2(NLEV-5:NLEV)
     !print*,AN(NLEV-5:NLEV)

     return
     end subroutine cmue_h


!EOC

!-----------------------------------------------------------------------
! Copyright by the GOTM-team under the GNU Public License - www.gnu.org
!-----------------------------------------------------------------------
