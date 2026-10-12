program hydrogen_reader
 use define
 use nlevel_cooling_mod, only: nlevel_cooling_add
 implicit none
 integer :: n
 call get_command_argument(1,par%atomic_dir)
 par%metal_excitation_model='invalid'
 n=nlevel_cooling_add('h',1)
 if(n<=0) error stop 'Hydrogen reader did not load the common atom'
 write(*,*) n
end program hydrogen_reader
