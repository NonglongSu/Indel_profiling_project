default:all
all:test1 test2

test1: Test_pseudo_seq.R/run_test.R
	Rscript --vanilla $<


test2: Test_update_gap.R/run_test.R
	Rscript --vanilla $<


.PHONY:test1 test2
