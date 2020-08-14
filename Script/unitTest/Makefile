default:all
all:    test1 test2

RSCRIPT = Rscript --vanilla 


###########################

test1: test_update_gap/run_test.R ../update_gap.R test_update_gap
	$(RSCRIPT) $< $(word 2,$^) $(word 3,$^)


test2: test_sw_gap/run_test.R ../sw_gap.R test_sw_gap
	$(RSCRIPT) $< $(word 2,$^) $(word 3,$^)



.PHONY:test1 test2
