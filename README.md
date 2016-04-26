# dns-pipe
Direct Numerical Simulation (DNS) of a pipe flow (CPL)

The present implementation allocates arrays of structures, 
which are allocated using several shared memory segments placed in the shared memory page.  The Linux kernel
managment of the shared memory page has predefined machine-dependent limits to the maximum size
of shared memory page kernel.shmall, to the maximum size of largest shared memory block kernel.shmmax and
to the maximum number of shared memory segments kernel.shmmni. 

In order to avoid apparent Out Of Memory errors, due to having reached one of those limits, set 
them to the maximum allowed value with 

```
echo "9999999999999" >/proc/sys/kernel/shmmax
echo "9999999999999" >/proc/sys/kernel/shmall
echo "9999999" >/proc/sys/kernel/shmmni
```

be careful! Those value may overflow, giving unexpected results! You can check your limit via

```
sysctl kernel.shm{max,all,mni}
```
