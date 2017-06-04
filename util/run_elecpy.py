import subprocess
import glob
import os

## Continue simulation for visualization

subprocess.call('python mahajan.py -i 400 -s //JALIFE/Recordings/SimulationResults/20170417-1')
subprocess.call('python elecpy.py -d //JALIFE/Recordings/SimulationResults/20170417-2 -m Mahajan')
#subprocess.call('python elecpy.py -r 0 -s E:/ExperimentData/20161118-1 -d //JALIFE/Recordings/SimulationResults/20170417-2 -m Mahajan')
#subprocess.call('python elecpy.py -p params/batch/param_-3.0_0.json -r 2540 -s E:/ExperimentData/20170111-30 -d E:/ExperimentData/20170111-30')
#subprocess.call('python elecpy.py -p params/batch/param_-3.0_1.json -r 2532 -s E:/ExperimentData/20170111-31 -d E:/ExperimentData/20170111-31')
#subprocess.call('python elecpy.py -p params/batch/param_-3.0_2.json -r 2573 -s E:/ExperimentData/20170111-32 -d E:/ExperimentData/20170111-32')

#subprocess.call('python elecpy.py -p params/batch/param_-2.0_0.json -r 2550 -s E:/ExperimentData/20170111-36 -d E:/ExperimentData/20170111-36')
#subprocess.call('python elecpy.py -p params/batch/param_-2.0_1.json -r 2542 -s E:/ExperimentData/20170111-37 -d E:/ExperimentData/20170111-37')
#subprocess.call('python elecpy.py -p params/batch/param_-2.0_2.json -r 2481 -s E:/ExperimentData/20170111-38 -d E:/ExperimentData/20170111-38')

#subprocess.call('python elecpy.py -p params/batch/param_-1.0_0.json -r 2557 -s E:/ExperimentData/20170111-42 -d E:/ExperimentData/20170111-42')
#subprocess.call('python elecpy.py -p params/batch/param_-1.0_1.json -r 2549 -s E:/ExperimentData/20170111-43 -d E:/ExperimentData/20170111-43')
#subprocess.call('python elecpy.py -p params/batch/param_-1.0_2.json -r 2488 -s E:/ExperimentData/20170111-44 -d E:/ExperimentData/20170111-44')

#subprocess.call('python elecpy.py -p params/batch/param_0.0_0.json -r 2577 -s E:/ExperimentData/20170111-48 -d E:/ExperimentData/20170111-48')
#subprocess.call('python elecpy.py -p params/batch/param_0.0_1.json -r 2569 -s E:/ExperimentData/20170111-49 -d E:/ExperimentData/20170111-49')
#subprocess.call('python elecpy.py -p params/batch/param_0.0_2.json -r 2505 -s E:/ExperimentData/20170111-50 -d E:/ExperimentData/20170111-50')

#subprocess.call('python elecpy.py -p params/batch/param_1.0_0.json -r 2513 -s E:/ExperimentData/20170111-15 -d E:/ExperimentData/20170111-15')
#subprocess.call('python elecpy.py -p params/batch/param_1.0_1.json -r 2505 -s E:/ExperimentData/20170111-16 -d E:/ExperimentData/20170111-16')
#subprocess.call('python elecpy.py -p params/batch/param_1.0_2.json -r 2540 -s E:/ExperimentData/20170111-17 -d E:/ExperimentData/20170111-17')

# subprocess.call('python elecpy.py -p params/batch/param_2.0_0.json -r 2526 -s E:/ExperimentData/20170111-21 -d E:/ExperimentData/20170111-21')
#subprocess.call('python elecpy.py -p params/batch/param_2.0_1.json -r 2518 -s E:/ExperimentData/20170111-22 -d E:/ExperimentData/20170111-22')
#subprocess.call('python elecpy.py -p params/batch/param_2.0_2.json -r 2557 -s E:/ExperimentData/20170111-23 -d E:/ExperimentData/20170111-23')

#subprocess.call('python elecpy.py -p params/batch/param_3.0_0.json -r 2537 -s E:/ExperimentData/20170111-27 -d E:/ExperimentData/20170111-27')
#subprocess.call('python elecpy.py -p params/batch/param_3.0_1.json -r 2529 -s E:/ExperimentData/20170111-28 -d E:/ExperimentData/20170111-28')
#subprocess.call('python elecpy.py -p params/batch/param_3.0_2.json -r 2570 -s E:/ExperimentData/20170111-29 -d E:/ExperimentData/20170111-29')


## Array stimulation

# subprocess.call('python elecpy.py -p params/batch/param_2400_22.json -r 2379 -s ./result/20170109-00 -d ./result/20170114-00')
# subprocess.call('python elecpy.py -p params/batch/param_2400_45.json -r 2379 -s ./result/20170109-00 -d ./result/20170114-01')
# subprocess.call('python elecpy.py -p params/batch/param_2400_90.json -r 2379 -s ./result/20170109-00 -d ./result/20170114-02')
# subprocess.call('python elecpy.py -p params/batch/param_2420_22.json -r 2399 -s ./result/20170109-00 -d ./result/20170114-03')
# subprocess.call('python elecpy.py -p params/batch/param_2420_45.json -r 2399 -s ./result/20170109-00 -d ./result/20170114-04')
# subprocess.call('python elecpy.py -p params/batch/param_2420_90.json -r 2399 -s ./result/20170109-00 -d ./result/20170114-05')
# subprocess.call('python elecpy.py -p params/batch/param_2440_22.json -r 2419 -s ./result/20170109-00 -d ./result/20170114-06')
# subprocess.call('python elecpy.py -p params/batch/param_2440_45.json -r 2419 -s ./result/20170109-00 -d ./result/20170114-07')
# subprocess.call('python elecpy.py -p params/batch/param_2440_90.json -r 2419 -s ./result/20170109-00 -d ./result/20170114-08')
# subprocess.call('python elecpy.py -p params/batch/param_2460_22.json -r 2439 -s ./result/20170109-00 -d ./result/20170114-09')
# subprocess.call('python elecpy.py -p params/batch/param_2460_45.json -r 2439 -s ./result/20170109-00 -d ./result/20170114-10')
# subprocess.call('python elecpy.py -p params/batch/param_2460_90.json -r 2439 -s ./result/20170109-00 -d ./result/20170114-11')


## S1-S2 stim

# subprocess.call('python elecpy.py -p params/temp.json -r 129 -s ./result/20170113-02 -d ./result/20170113-02')


## single point stim

# subprocess.call('python elecpy.py -p params/temp.json -r 0 -s ./result/20170113-00 -d ./result/20170113-01')


## Sequential stimulation on a line

# subprocess.call('python elecpy.py -p params/batch/param_1.0_0.json -r 2406 -s ./result/20170109-00 -d ./result/20170111-15')
# subprocess.call('python elecpy.py -p params/batch/param_1.0_1.json -r 2400 -s ./result/20170109-00 -d ./result/20170111-16')
# subprocess.call('python elecpy.py -p params/batch/param_1.0_2.json -r 2440 -s ./result/20170109-00 -d ./result/20170111-17')

# subprocess.call('python elecpy.py -p params/batch/param_1.5_0.json -r 2413 -s ./result/20170109-00 -d ./result/20170111-18')
# subprocess.call('python elecpy.py -p params/batch/param_1.5_1.json -r 2403 -s ./result/20170109-00 -d ./result/20170111-19')
# subprocess.call('python elecpy.py -p params/batch/param_1.5_2.json -r 2450 -s ./result/20170109-00 -d ./result/20170111-20')

# subprocess.call('python elecpy.py -p params/batch/param_2.0_0.json -r 2418 -s ./result/20170109-00 -d ./result/20170111-21')
# subprocess.call('python elecpy.py -p params/batch/param_2.0_1.json -r 2408 -s ./result/20170109-00 -d ./result/20170111-22')
# subprocess.call('python elecpy.py -p params/batch/param_2.0_2.json -r 2458 -s ./result/20170109-00 -d ./result/20170111-23')

# subprocess.call('python elecpy.py -p params/batch/param_2.5_0.json -r 2423 -s ./result/20170109-00 -d ./result/20170111-24')
# subprocess.call('python elecpy.py -p params/batch/param_2.5_1.json -r 2413 -s ./result/20170109-00 -d ./result/20170111-25')
# subprocess.call('python elecpy.py -p params/batch/param_2.5_2.json -r 2465 -s ./result/20170109-00 -d ./result/20170111-26')

# subprocess.call('python elecpy.py -p params/batch/param_3.0_0.json -r 2428 -s ./result/20170109-00 -d ./result/20170111-27')
# subprocess.call('python elecpy.py -p params/batch/param_3.0_1.json -r 2419 -s ./result/20170109-00 -d ./result/20170111-28')
# subprocess.call('python elecpy.py -p params/batch/param_3.0_2.json -r 2472 -s ./result/20170109-00 -d ./result/20170111-29')

# subprocess.call('python elecpy.py -p params/batch/param_-3.0_0.json -r 2431 -s ./result/20170109-00 -d ./result/20170111-30')
# subprocess.call('python elecpy.py -p params/batch/param_-3.0_1.json -r 2421 -s ./result/20170109-00 -d ./result/20170111-31')
# subprocess.call('python elecpy.py -p params/batch/param_-3.0_2.json -r 2473 -s ./result/20170109-00 -d ./result/20170111-32')

# subprocess.call('python elecpy.py -p params/batch/param_-2.5_0.json -r 2437 -s ./result/20170109-00 -d ./result/20170111-33')
# subprocess.call('python elecpy.py -p params/batch/param_-2.5_1.json -r 2426 -s ./result/20170109-00 -d ./result/20170111-34')
# subprocess.call('python elecpy.py -p params/batch/param_-2.5_2.json -r 2477 -s ./result/20170109-00 -d ./result/20170111-35')

# subprocess.call('python elecpy.py -p params/batch/param_-2.0_0.json -r 2442 -s ./result/20170109-00 -d ./result/20170111-36')
# subprocess.call('python elecpy.py -p params/batch/param_-2.0_1.json -r 2430 -s ./result/20170109-00 -d ./result/20170111-37')
# subprocess.call('python elecpy.py -p params/batch/param_-2.0_2.json -r 2381 -s ./result/20170109-00 -d ./result/20170111-38')

# subprocess.call('python elecpy.py -p params/batch/param_-1.5_0.json -r 2445 -s ./result/20170109-00 -d ./result/20170111-39')
# subprocess.call('python elecpy.py -p params/batch/param_-1.5_1.json -r 2432 -s ./result/20170109-00 -d ./result/20170111-40')
# subprocess.call('python elecpy.py -p params/batch/param_-1.5_2.json -r 2384 -s ./result/20170109-00 -d ./result/20170111-41')

# subprocess.call('python elecpy.py -p params/batch/param_-1.0_0.json -r 2450 -s ./result/20170109-00 -d ./result/20170111-42')
# subprocess.call('python elecpy.py -p params/batch/param_-1.0_1.json -r 2439 -s ./result/20170109-00 -d ./result/20170111-43')
# subprocess.call('python elecpy.py -p params/batch/param_-1.0_2.json -r 2388 -s ./result/20170109-00 -d ./result/20170111-44')

# subprocess.call('python elecpy.py -p params/batch/param_-0.5_0.json -r 2459 -s ./result/20170109-00 -d ./result/20170111-45')
# subprocess.call('python elecpy.py -p params/batch/param_-0.5_1.json -r 2447 -s ./result/20170109-00 -d ./result/20170111-46')
# subprocess.call('python elecpy.py -p params/batch/param_-0.5_2.json -r 2398 -s ./result/20170109-00 -d ./result/20170111-47')

# subprocess.call('python elecpy.py -p params/batch/param_0.0_0.json -r 2379 -s ./result/20170109-00 -d ./result/20170111-48')
# subprocess.call('python elecpy.py -p params/batch/param_0.0_1.json -r 2462 -s ./result/20170109-00 -d ./result/20170111-49')
# subprocess.call('python elecpy.py -p params/batch/param_0.0_2.json -r 2405 -s ./result/20170109-00 -d ./result/20170111-50')

# subprocess.call('python elecpy.py -p params/batch/param_0.5_0.json -r 2393 -s ./result/20170109-00 -d ./result/20170111-51')
# subprocess.call('python elecpy.py -p params/batch/param_0.5_1.json -r 2382 -s ./result/20170109-00 -d ./result/20170111-52')
# subprocess.call('python elecpy.py -p params/batch/param_0.5_2.json -r 2422 -s ./result/20170109-00 -d ./result/20170111-53')

## Simultaneous stimulation on a line

# files = glob.glob('./params/batch/*.json')
# for i, f in enumerate(files):
# 	cmd = 'python elecpy.py -p {0} -r 2350 -s ./result/20170109-00 -d ./result/20170111-{1:0>2}'.format(f, i+1)
# 	print cmd
# 	os.system(cmd)

## Sequential stimulation on a line
# subprocess.call('python elecpy.py -p params/temp.json -r 2414 -s ./result/temp -d ./result/20170110-00')
# subprocess.call('python elecpy.py -p params/temp.json -r 2414 -s ./result/20170109-00 -d ./result/temp')

# subprocess.call('python elecpy.py -p params/temp.json -r 2009 -s ./result/20161124-1 -d ./result/temp')

# subprocess.call('python elecpy.py -p params/temp.json -r 2409 -s ./result/20170109-00 -d ./result/temp')

# subprocess.call('python elecpy.py -p params/2420/25_0.json -r 2419 -s ./result/20170109-00 -d ./result/20170109-01')

## Re-start from pinning state
# subprocess.call('python elecpy.py -p params/temp.json -r 2019 -s ./result/20161124-1 -d ./result/20170109-00')

## 150x200 trial
# subprocess.call('python elecpy.py -p params/temp.json -r 0 -s ./result/20170107-01 -d ./result/20170108-01')
# subprocess.call('python elecpy.py -p params/temp.json -d ./result/temp')

# subprocess.call('python elecpy.py -p params/2020/75_0.json -r 2019 -s ./result/1124-1 -d ./result/1124-3')
# subprocess.call('python elecpy.py -p params/2020/75_1.json -r 2019 -s ./result/1124-1 -d ./result/1124-4')
# subprocess.call('python elecpy.py -p params/2020/75_2.json -r 2019 -s ./result/1124-1 -d ./result/1124-5')
# subprocess.call('python elecpy.py -p params/2020/75_3.json -r 2019 -s ./result/1124-1 -d ./result/1124-6')
# subprocess.call('python elecpy.py -p params/2020/75_4.json -r 2019 -s ./result/1124-1 -d ./result/1124-7')
# subprocess.call('python elecpy.py -p params/2020/75_5.json -r 2019 -s ./result/1124-1 -d ./result/1124-8')

# subprocess.call('python elecpy.py -p params/2040/75_0.json -r 2039 -s ./result/1124-1 -d ./result/1124-9')
# subprocess.call('python elecpy.py -p params/2040/75_1.json -r 2039 -s ./result/1124-1 -d ./result/1124-10')
# subprocess.call('python elecpy.py -p params/2040/75_2.json -r 2039 -s ./result/1124-1 -d ./result/1124-11')
# subprocess.call('python elecpy.py -p params/2040/75_3.json -r 2039 -s ./result/1124-1 -d ./result/1124-12')
# subprocess.call('python elecpy.py -p params/2040/75_4.json -r 2039 -s ./result/1124-1 -d ./result/1124-13')
# subprocess.call('python elecpy.py -p params/2040/75_5.json -r 2039 -s ./result/1124-1 -d ./result/1124-14')

# subprocess.call('python elecpy.py -p params/2060/75_0.json -r 2059 -s ./result/1124-1 -d ./result/1124-15')
# subprocess.call('python elecpy.py -p params/2060/75_1.json -r 2059 -s ./result/1124-1 -d ./result/1124-16')
# subprocess.call('python elecpy.py -p params/2060/75_2.json -r 2059 -s ./result/1124-1 -d ./result/1124-17')
# subprocess.call('python elecpy.py -p params/2060/75_3.json -r 2059 -s ./result/1124-1 -d ./result/1124-18')
# subprocess.call('python elecpy.py -p params/2060/75_4.json -r 2059 -s ./result/1124-1 -d ./result/1124-19')
# subprocess.call('python elecpy.py -p params/2060/75_5.json -r 2059 -s ./result/1124-1 -d ./result/1124-20')

# subprocess.call('python elecpy.py -p params/2080/75_0.json -r 2079 -s ./result/1124-1 -d ./result/1124-21')
# subprocess.call('python elecpy.py -p params/2080/75_1.json -r 2079 -s ./result/1124-1 -d ./result/1124-22')
# subprocess.call('python elecpy.py -p params/2080/75_2.json -r 2079 -s ./result/1124-1 -d ./result/1124-23')
# subprocess.call('python elecpy.py -p params/2080/75_3.json -r 2079 -s ./result/1124-1 -d ./result/1124-24')
# subprocess.call('python elecpy.py -p params/2080/75_4.json -r 2079 -s ./result/1124-1 -d ./result/1124-25')
# subprocess.call('python elecpy.py -p params/2080/75_5.json -r 2079 -s ./result/1124-1 -d ./result/1124-26')

# subprocess.call('python elecpy.py -p params/2100/75_0.json -r 2099 -s ./result/1124-1 -d ./result/1124-27')
# subprocess.call('python elecpy.py -p params/2100/75_1.json -r 2099 -s ./result/1124-1 -d ./result/1124-28')
# subprocess.call('python elecpy.py -p params/2100/75_2.json -r 2099 -s ./result/1124-1 -d ./result/1124-29')
# subprocess.call('python elecpy.py -p params/2100/75_3.json -r 2099 -s ./result/1124-1 -d ./result/1124-30')
# subprocess.call('python elecpy.py -p params/2100/75_4.json -r 2099 -s ./result/1124-1 -d ./result/1124-31')
# subprocess.call('python elecpy.py -p params/2100/75_5.json -r 2099 -s ./result/1124-1 -d ./result/1124-32')
