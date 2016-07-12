#include"cppdefs.h"
!-----------------------------------------------------------------------
!BOP
!
! !ROUTINE: Update dimensionless alpha's\label{sec:alpha}
!
! !INTERFACE:
   subroutine alpha_mnb(nlev,NN,SS,SSU,SSV,SSS,SSUS,SSVS)
!
! !DESCRIPTION:
! This subroutine updates the dimensionless numbers $\alpha_M$, $\alpha_N$,
! and $\alpha_b$ according to \eq{alphaMN}. Note that according to \eq{Nbar}
! and \eq{NbarVertical} the following identities are valid
! \begin{equation}
!  \label{alphaIdentities}
!    \alpha_M = \overline{S}^2 \comma
!    \alpha_N = \overline{N}^2 \comma
!    \alpha_b = \overline{T}   \point
! \end{equation}
!
!
! !USES:
  use turbulence,  only:     tke,eps,kb,l
  use turbulence,  only:     as,an,at,astk,av
  IMPLICIT NONE
!
! !INPUT PARAMETERS:
  integer,  intent(in)      :: nlev
  REALTYPE, intent(in)      :: NN(0:nlev),SS(0:nlev), SSS(0:nlev)
  REALTYPE, intent(in)      :: SSU(0:nlev), SSV(0:nlev)
  REALTYPE, intent(in)      :: SSUS(0:nlev), SSVS(0:nlev)

!
! !REVISION HISTORY:
!  Original author(s): Lars Umlauf
!
!EOP
!-----------------------------------------------------------------------
! !LOCAL VARIABLES:
  integer              :: i
  REALTYPE             :: tau2

!-----------------------------------------------------------------------
!BOC

  do i=0,nlev
     tau2   = tke(i)*tke(i) / ( eps(i)*eps(i) )
     as(i)  = tau2 * SS(i)
     an(i)  = tau2 * NN(i)
     astk(i)  = tau2 * SSS(i)
     av(i)  = tau2 * sqrt(SSU(i)*SSUS(i)+SSV(i)*SSVS(i))
     at(i)  = tke(i)/eps(i) * kb(i)/eps(i)

!    clip negative values
     as(i) = max(as(i),1.e-10*_ONE_)
     at(i) = max(at(i),1.e-10*_ONE_)
     astk(i) = max(astk(i),_ZERO_)
     av(i) = max(av(i),1.e-10*_ZERO_)
  end do
  
  return
end subroutine alpha_mnb

!EOC

!-----------------------------------------------------------------------
! Copyright by the GOTM-team under the GNU Public License - www.gnu.org
!-----------------------------------------------------------------------
