import os
import sys
import random

if __name__ == '__main__':
    if len(sys.argv) < 4:
        exit(0)
    
    board_size = int(sys.argv[1])
    target_gen = int(sys.argv[2])
    ghost_size = int(sys.argv[3])
    num_iter = 25

    sample_dir = 'my_sample'
    os.makedirs(sample_dir, exist_ok=True)
    
    for i in range(num_iter):
        input_name  = f'input{i}.txt'
        output_name = f'output{i}.txt'
        
        input_file = open(f'{sample_dir}/{input_name}', 'w')
        input_file.write(f'{board_size}\n{target_gen}\n{ghost_size}\n')
        for p in range(board_size):
            for q in range(board_size):
                c = random.randint(0, 12) % 4
                if c == 0:
                    input_file.write('#')
                else:
                    input_file.write('.')
            input_file.write('\n')
        input_file.close()

        os.system(f'./conway_single < ./{sample_dir}/{input_name} > ./{sample_dir}/{output_name}')