program main
use feature_main
type(feature_manager) :: All_result
INTEGER I,leibie_temp,Ti,Tj,Tk,Numi
INTEGER Size_m,Size_n,Size_l
leibie_temp=0

!�����ļ�: files,�����������һ�У����ļ�����������ÿ��POSCAR��Ӧ������

call ALL_result%max_atoms()

call ALL_result%check_level(All_result%POSCAR_name(All_result%num_of_max_atoms_file),leibie_temp)
ALL_result%num_of_lattice(3)=leibie_temp
!print *,leibie_temp
call ALL_result%check_level_m(All_result%POSCAR_name(All_result%num_of_max_atoms_file),leibie_temp)
ALL_result%num_of_lattice(1)=leibie_temp
call ALL_result%check_level_n(All_result%POSCAR_name(All_result%num_of_max_atoms_file),leibie_temp)
ALL_result%num_of_lattice(2)=leibie_temp


ALLOCATE(ALL_result%Feature(ALL_result%num_of_lattice(1),ALL_result%num_of_lattice(2)&
,ALL_result%num_of_lattice(3)+1,ALL_result%num_of_file))
!������Se��Sԭ�ӣ���˱�MoҪ����һ��
 
call ALL_result%read_poscar()
Size_m=size(ALL_result%Feature,1)
Size_n=size(ALL_result%Feature,2)
Size_l=size(ALL_result%Feature,3)

!������
open(unit=1098,file='Features.dat')
DO Numi=1,All_result%num_of_file
    write(1098,*) ALL_result%Feature(:,:,:,Numi)
ENDDO
close(1098)
!!only for check
!I=1
!DO Ti=1,2
!    DO Tj=1,4
!        DO Tk=1,2
!         ALL_result%Feature(Ti,Tj,Tk,1)=I   
!       I=I+1
!        ENDDO
!    ENDDO
!ENDDO
!
!
!open(unit=1098,file='FeaCheck.dat')
!    write(1098,*) ALL_result%Feature(:,:,:,1)
!close(1098)
end program main