 subroutine mark_layers(mark_of_centers,mark_of_surrounds,&!now mark_of_surrounds are true for all atoms of one sort or all sorts
 &z1_for_centres,z2_for_centres,sort_of_centers,sort_of_surronds,vec_orht_to_layers,&
 &r_at,n_at,i_sort_at,sizex,sizey,sizez)
 use defs_basis
 include 'max_at.h'
 implicit none

 logical(lgt),intent(out) :: mark_of_centers(max_at),mark_of_surrounds(max_at)
 real(dp),intent(in) :: r_at(3,max_at),n_at,sizex,sizey,sizez,z1_for_centres,z2_for_centres
 integer,intent(in) :: i_sort_at(max_at),sort_of_centers,sort_of_surronds,vec_orht_to_layers(3)


  integer :: i
 real(dp) :: zmin,z1,z2

 mark_of_centers=.false.
 mark_of_surrounds=.false.
 zmin=10000.
 do i=1,n_at
 if(r_at(3,i)<zmin) zmin=r_at(3,i)
 enddo
 z1=zmin+z1_for_centres; z2=zmin+z2_for_centres


 do i=1,n_at
 if(r_at(3,i)>z1.and.r_at(3,i)<z2)then
 if(i_sort_at(i)==sort_of_centers.or.sort_of_centers==0)) mark_of_centers(i)=.true.
 endif
 enddo

 end subroutine mark_layers