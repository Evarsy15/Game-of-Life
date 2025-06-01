import os

exe_list = ['conway_single']
num_proc = 4

if __name__ == '__main__':
    for exe in exe_list:
        if exe == 'conway_single':
            for i in range(1, 6):
                os.system(f'./{exe} < ./sample/input{i}.txt > my_output{i}.txt')
                os.system(f'diff -bwi my_output{i}.txt ./sample/output{i}.txt > my_log{i}.txt')
        elif exe in ['conway_2d_mesh', 'conway_1d_mesh']:
            for i in range(1, 6):
                os.system(f'mpirun -np {num_proc} ./{exe} < ./sample/input{i}.txt > my_output{i}.txt')
                os.system(f'diff -bwi my_output{i}.txt ./sample/output{i}.txt > my_log{i}.txt')