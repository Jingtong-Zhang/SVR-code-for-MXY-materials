module feature_main
    

implicit none
private
type,public :: feature_manager

character(len=512) File_name 
character(len=512),allocatable :: POSCAR_name(:)
integer num_of_file !�ܹ��ж��ٸ�POSCAR
INTEGER num_of_max_atoms !���ԭ����
INTEGER num_of_max_Mo
INTEGER num_of_max_atoms_file !�������ԭ�������ļ��ţ���Ӧ��POSCAR_name(num_of_max_atoms_file)��
INTEGER num_of_lattice(3) !N*M*L�Ľṹ
REAL,allocatable :: Feature(:,:,:,:) !���������feature��������ʲôԭ��
contains
procedure :: read_poscar
procedure :: max_atoms !�������ԭ�������ӳ���
procedure :: check_level !�ж���N*M*L�Ľṹ���õ�L��ֵ
procedure :: check_level_m
procedure :: check_level_n
procedure :: paixu   !ð�����򷨰�position��С��������

end type feature_manager
contains
subroutine read_poscar(self)
!�������е�POSCAR�ļ����������ǵĽ��������self���type����ȥ
class(feature_manager),intent(inout) :: self
REAL,allocatable :: location(:,:,:)
REAL,allocatable :: position_l(:,:)
character(len=512) cFile_test
INTEGER total_atom
INTEGER I,J
INTEGER Ti,Tj,Tk,Uei,Uej,Uek,Wri
character(len=512) DorS
INTEGER Size_m,Size_n,Size_l !�������󾧰���size
INTEGER judge_m,judge_n,judge_l!����ǵ�ǰ������size
INTEGER beilv_m,beilv_n,beilv_l
INTEGER atoms(3)
Size_m=size(self%Feature,1)
Size_n=size(self%Feature,2)
Size_l=size(self%Feature,3)

ALLOCATE(location(Size_m,Size_n,Size_l))
!position_l����洢ԭ�ӵ������Լ���ԭ�ӵ�����

DO I=1,self%num_of_file
call self%check_level_m(self%POSCAR_name(I),judge_m)
call self%check_level_n(self%POSCAR_name(I),judge_n)
call self%check_level(self%POSCAR_name(I),judge_l)
judge_l=judge_l+1
beilv_m=Size_m/judge_m
beilv_n=Size_n/judge_n
beilv_l=Size_l/judge_l
!print *,self%POSCAR_name(I)
open(unit=107,file=self%POSCAR_name(I))
read(107,*)
read(107,*)
read(107,*)
read(107,*)
read(107,*)
read(107,*)
read(107,*) atoms 
!
total_atom=atoms(2)+atoms(3)
if (ALLOCATED(position_l))then
DEALLOCATE(position_l)
!����Ѿ������ˣ��ǾͰ����ͷŵ�
end if 
ALLOCATE(position_l(4,total_atom))

read(107,*) DorS

DorS=Trim(adjustL(DorS))
if (DorS=='Direct')then
else
read(107,*) !Selectritive dynamics
end if 
DO J=1,atoms(1)
    read(107,*) !�������е�Moԭ��
ENDDO
DO J=1,atoms(2)
read(107,*) position_l(1:3,J)
position_l(4,J)=1 !1����Seԭ��
ENDDO
!print *,atoms(2),atoms(3)
DO J=(atoms(2)+1),(atoms(2)+atoms(3))
   ! print *,J
read(107,*) position_l(1:3,J)
position_l(4,J)=-1 !1����Sԭ��
ENDDO
close(107)
!����
CALL self%paixu(position_l,total_atom)
!����֮��position_l����xyz��˳���С��������

!д�뵽feature����ȥ
Uei=1
Uej=1
Uek=1
write(cFile_test,*) I
open(unit=1093,file=''//Trim(adjustL(cFile_test))//'')
DO Ti=1,Size_m
    DO Tj=1,Size_n
        DO Tk=1,Size_l
            Wri=judge_l*judge_n*(Uei-1)+judge_l*(Uej-1)+Uek
            self%Feature(Ti,Tj,Tk,I)=position_l(4,Wri)
            write(1093,*) position_l(1:4,Wri),Uei,Uej,Uek,Wri
          !  if(I==3)then
           ! print *,position_l(4,Wri)
            !end if 
            Uek=Uek+1
            if(Uek>judge_l)then
                Uek=1
                if (Tk==Size_l)then
                Uej=Uej+1
                if(Uej>judge_n)then
                Uej=1
                if(Tj==Size_n)then
                Uei=Uei+1
                if (Uei>judge_m)then
                Uei=1
                end if 
                end if 
                end if 
                end if 
            end if 
            
        ENDDO
    ENDDO
ENDDO
!��position_l�������,���������û����


close(1093)
ENDDO





end subroutine

subroutine paixu(self,Position_all,siza)
class(feature_manager),intent(inout) :: self
INTEGER siza
REAL Position_all(4,siza)
REAL Temp1(4),Temp2(4)
INTEGER I,J,K
REAL Finish

!�ȶ����е�x��������
Finish=1
do while(Finish>0) 
Finish=-1
DO J=1,siza-1
Temp1=Position_all(:,J)
Temp2=Position_all(:,J+1)
if (Temp1(1)>Temp2(1))then
Position_all(:,J+1)=Temp1
Position_all(:,J)=Temp2
Finish=1
end if 

ENDDO
enddo 
!����y��������
Finish=1
do while(Finish>0) 
Finish=-1
DO J=1,siza-1
Temp1=Position_all(:,J)
Temp2=Position_all(:,J+1)
if ((Temp1(2)>Temp2(2)).and.(Temp1(1)>=Temp2(1)))then
Position_all(:,J+1)=Temp1
Position_all(:,J)=Temp2
Finish=1
end if 

ENDDO
enddo 
!����z��������
Finish=1
do while(Finish>0) 
Finish=-1
DO J=1,siza-1
Temp1=Position_all(:,J)
Temp2=Position_all(:,J+1)
if ((Temp1(1)>=Temp2(1)).and.(Temp1(2)>=Temp2(2)).and.(Temp1(3)>Temp2(3)))then
Position_all(:,J+1)=Temp1
Position_all(:,J)=Temp2
Finish=1
end if 

ENDDO
enddo 


end subroutine paixu 

subroutine max_atoms(self)
class(feature_manager),intent(inout) :: self
INTEGER max_atom_I,atom_temp(4)
character(len=512) cFile_max_atoms
open(unit=101,file='files')
read(101,*) self%num_of_file
ALLOCATE(self%POSCAR_name(self%num_of_file))
!����files���������������ļ����������POSCAR�ļ���

atom_temp=0
self%num_of_max_atoms=0
self%num_of_max_atoms_file=0
self%num_of_max_Mo=0
DO max_atom_I=1,self%num_of_file
read(101,*) cFile_max_atoms  !����POSCAR������
self%POSCAR_name(max_atom_I)=Trim(adjustL(cFile_max_atoms))
!�Ƴ������еĿո񣬲��ҽ���д�뵽POSCAR_name����ȥ
open(unit=102,file=self%POSCAR_name(max_atom_I))
read(102,*)
read(102,*)
read(102,*)
read(102,*)
read(102,*)
read(102,*)
read(102,*) atom_temp(1:3)
atom_temp(4)=sum(atom_temp(1:3))

if(atom_temp(4)>self%num_of_max_atoms)then
self%num_of_max_atoms=atom_temp(4)
self%num_of_max_atoms_file=max_atom_I
self%num_of_max_Mo=atom_temp(1)
end if 
close(102)
ENDDO

close(101)

end subroutine max_atoms

subroutine check_level(self,cfile_l,judge2)
class(feature_manager),intent(inout) :: self
REAL Temp_acell(3),Hight,Temp1,Temp2
REAL,allocatable :: location(:)
character(len=512) cfile_l
INTEGER Temp_int(3),Temp_I
INTEGER judge2
REAL Hightest,Lowest,Diff
character(len=512) DorS
!�ж�L
cfile_l=Trim(adjustL(cfile_l))
open(unit=103,file=cfile_l)
!�򿪾������ԭ�������ļ�
read(103,*)
read(103,*)
read(103,*)
read(103,*)
read(103,*) Temp_acell
Hight=Temp_acell(3)
read(103,*)
read(103,*) Temp_int
ALLOCATE(location(Temp_int(1))) 
read(103,*) DorS
DorS=Trim(adjustL(DorS))
if (DorS=='Direct')then
else
read(103,*) !Selectritive dynamics
end if 

DO Temp_I=1,Temp_int(1)!�������е�Moԭ��
read(103,*) Temp1,Temp2,location(Temp_I)
location(Temp_I)=location(Temp_I)*Hight
ENDDO
Hightest=maxval(location)
Lowest=minval(location)
Diff=Hightest-Lowest
if (Diff>1.0)then
print *,'More than one layer is found, the results is untrustable'
else
judge2=1
end if 
close(103)
end subroutine

subroutine check_level_m(self,file_m,judge2)
class(feature_manager),intent(inout) :: self
character(len=512) DorS
character(len=512) file_m
INTEGER judge2
INTEGER I,J
INTEGER Leibie,ATOMS(3)
REAL,allocatable :: location(:)
REAL Temp(3),Dif
INTEGER Judge
ALLOCATE(location(self%num_of_max_Mo)) 
location=2
file_m=Trim(adjustL(file_m))
open(unit=105,file=file_m)
DO I=1,6
    read(105,*)
ENDDO
read(105,*) ATOMS
read(105,*) DorS
DorS=Trim(adjustL(DorS))
if (DorS=='Direct')then
else
read(105,*) !Selectritive dynamics
end if 
Leibie=0
read(105,*) Temp
location(1)=Temp(1)
Leibie=Leibie+1 
DO I=2,ATOMS(1)
!print *,I
read(105,*) Temp
!�������ж�Temp(1)��location����ÿһ����ľ���
!�����location���������һ��ֵС��ĳ��ֵ����Leibie+1
Judge=0
DO J=1,ATOMS(1)
Dif=abs(Temp(1)-location(J))
if (Dif<(1/(ATOMS(1)*1.0)))then
Judge=1
end if 
ENDDO
if (Judge==1)then
Leibie=Leibie+1
location(Leibie)=Temp(1)
end if 
ENDDO

close(105)

judge2=Leibie

end subroutine check_level_m
subroutine check_level_n(self,file_n,judge2)
class(feature_manager),intent(inout) :: self
character(len=512) DorS
character(len=512) file_n
INTEGER judge2
INTEGER I,J
INTEGER ATOMS(3)
INTEGER Leibie
REAL,allocatable :: location(:)
REAL Temp(3),Dif
INTEGER Judge
ALLOCATE(location(self%num_of_max_Mo)) 
location=2
file_n=Trim(adjustL(file_n))
open(unit=105,file=file_n)
DO I=1,6
    read(105,*)
ENDDO
read(105,*) ATOMS
read(105,*) DorS
DorS=Trim(adjustL(DorS))
if (DorS=='Direct')then
else
read(105,*) !Selectritive dynamics
end if 
Leibie=0
read(105,*) Temp
location(1)=Temp(2)
Leibie=Leibie+1 
DO I=2,ATOMS(1)
read(105,*) Temp
!�������ж�Temp(1)��location����ÿһ����ľ���
!�����location���������һ��ֵС��ĳ��ֵ����Leibie+1
Judge=0
DO J=1,ATOMS(1)
Dif=abs(Temp(2)-location(J))
if (Dif<(1/(ATOMS(1)*1.0)))then
Judge=1
end if 
ENDDO
if (Judge==1)then
Leibie=Leibie+1
location(Leibie)=Temp(2)
end if 
ENDDO

close(105)
judge2=Leibie
end subroutine check_level_n





end module feature_main

