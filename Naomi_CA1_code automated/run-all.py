#!/usr/bin/python3

import os
from colorama import Fore

def compile_and_run(program_path):
	"""
	Compiles and runs an instance of the 'Naomi_CA_1' program located in the 
	'program_path' directory. 
	"""

	os.chdir(program_path)
	print(Fore.GREEN + f'[INFO] Compiling {program_path}' + Fore.RESET)
	os.system(f'make -f Makefile')
	print(Fore.GREEN + f'[INFO] Running {program_path}' + Fore.RESET)
	os.system(f'chmod +x mySim')
	os.system(f'./mySim > results/terminal_output.txt')
	os.chdir('..')


def process_results(results_path, results_dir):
	"""
	Copies results folder from 'results_path' into the results directory located
	in 'results_dir'.
	"""
	print(Fore.GREEN + f'[INFO] Transferring results from {results_path}' +\
		Fore.RESET)
	os.system(f'cp -r {results_path} ../../{results_dir}')

def process_dir(folder_path, processed_path, results_path):
	"""
	Processes all program folders with a name starting with 'prc_', and copies
	the results into 'results_path'.
	"""

	print(Fore.GREEN + f'[INFO] Processing \'{folder_path}\', storing results in \'{processed_path}\'' + Fore.RESET)
	os.chdir(folder_path)
	programs = os.listdir()
	for program in programs:
		if program[0] == '.':
			continue
		compile_and_run(program)
		process_results(f'\'{program}\'/\'{results_path}\'', results_path)
		print(Fore.GREEN + f'[INFO] Moving {program} to finished folder' + Fore.RESET)
		os.system(f'mv \'{program}\' ../../{processed_path}')
		print(Fore.GREEN + f'[INFO] Renaming ../../{results_path}/results to ../../{results_path}/\'{program} results\'')
		os.system(f'mv ../../{results_path}/results ../../{results_path}/\'{program} results\'')
	print(Fore.GREEN + f'[INFO] Finished proceesing \'{folder_path}\', results located in directory \'{results_path}\'' + Fore.RESET)
	os.chdir('..')

if __name__ == '__main__':
	process_dir('src/process/', 'src/finished/', 'results')