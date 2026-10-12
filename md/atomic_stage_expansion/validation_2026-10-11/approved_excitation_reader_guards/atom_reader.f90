program atom_reader
 use nlevel_mod, only: nlevel_atom_type,nlevel_load
 implicit none
 type(nlevel_atom_type)::a
 character(len=1024)::file
 logical::ok
 call get_command_argument(1,file)
 call nlevel_load(trim(file),a,ok)
 if(.not.ok) error stop 'Missing atom'
 write(*,*)a%nlev
end program atom_reader
