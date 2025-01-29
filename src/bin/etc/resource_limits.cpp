//+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-
// Copyright (C) 2024 Edward E. Meyer
// This file is part of Flaps; Flaps is free software: you can redistribute
// it and/or modify it under the terms of the GNU General Public License.
// See the file COPYING in the root directory.
// Flaps is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
//+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-

// Print resource limits that can be changed with setrlimit
// Ref: W. Richard Stevens, "Advanced Programming in the Unix Environment",
//      Addison-Wesley, 1992

#include	<cstdio>
#include	<sys/resource.h>
#include	<cstdlib>

static void	pr_limits(char const*, int);

int
main(void) {
	printf(" resource         current     maximum\n");
	pr_limits("core file",RLIMIT_CORE);
	pr_limits("cpu time", RLIMIT_CPU);
	pr_limits("data segment", RLIMIT_DATA);
	pr_limits("file size", RLIMIT_FSIZE);
#ifdef	RLIMIT_MEMLOCK
	pr_limits("locked-in-memory", RLIMIT_MEMLOCK);
#endif
#ifdef	RLIMIT_NOFILE	/* SVR4 name */
	pr_limits("open files",RLIMIT_NOFILE);
#endif
#ifdef	RLIMIT_OFILE	/* 44BSD name */
	pr_limits("open files(bsd)",RLIMIT_OFILE);
#endif
#ifdef	RLIMIT_NPROC
	pr_limits("child procs", RLIMIT_NPROC);
#endif
#ifdef	RLIMIT_RSS
	pr_limits("resident set size (RSS)", RLIMIT_RSS);
#endif
	pr_limits("stack size", RLIMIT_STACK);
#ifdef	RLIMIT_VMEM
	pr_limits("mmap size", RLIMIT_VMEM);
#endif
	exit(0);
}

static void
pr_limits(char const* name, int resource)
{
	struct rlimit	limit;

	if (getrlimit(resource, &limit) < 0)
		fprintf(stderr, "getrlimit error for %s", name);
	printf("%-23s  ", name);
	if (limit.rlim_cur == RLIM_INFINITY)
		printf("(infinite)  ");
	else
		printf("%10ld  ", limit.rlim_cur);
	if (limit.rlim_max == RLIM_INFINITY)
		printf("(infinite)\n");
	else
		printf("%10ld\n", limit.rlim_max);
}
