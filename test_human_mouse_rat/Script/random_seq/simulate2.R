library(seqinr)
library(stringr)

# fix issues with seqinr
source("write_fasta.R")

# maximum length of a M before glue starts
G = 99

gaps = scan("gaps.txt",character(0),quiet=TRUE)
cigar = gaps %>% str_extract_all(pattern="\\d+\\w") %>% lapply(., function(x) {
    lens = as.integer(x %>% str_extract_all(pattern="^\\d+"))
    chars = x %>% str_extract_all(pattern="\\w$")
    a = setNames(lens,chars)
    a
})

p_nucs = c("A"=0.308,"C"=0.185,"G"=0.199,"T"=0.308)
nucs = names(p_nucs)

# input = "random_seq/sample.fa"

simulate_two_main = function(input,ouDir) {
    seqs = read.fasta(input, set.attributes=FALSE, forceDNAtolower=FALSE)
    repeat {
        GG = G
        sims = seqs
        A = sims$mouse
        B = sims$rat
        stopifnot(length(A) == length(B), length(A) %% 3 == 0)

        # sample from cigar strings
        g = sample(cigar, 1)[[1]]
        # Identify any flexible width matches
        bDM = names(g) %in% c("M","D")
        repeat {
            bF = names(g) == "M" & g >= GG
            if(sum(bF) != 0) {
                break
            }
            GG = (GG %/% 6)*3
        }
        # save phase
        p = g %% 3
        # total committed length
        w = sum(c(p[bF]+GG,g[bDM & !bF]))
        if(w > length(A)) {
            # this cigar is too big for our sequence, so try again
            next
        }
        # number of glue sections
        n = sum(bF)
        # sample glue widths from the remaining length
        stopifnot((length(A)-w) %% 3 == 0)
        r = (length(A)-w) %/% 3
        s = round(r*sort(runif(n-1)))
        x = c(s,r)-c(0,s)
        # update widths
        g[bF] = GG+p[bF]+3*x
        stopifnot(g %% 3 == p)
        stopifnot(sum(g[bDM]) == length(A))

        # update alignments
        pos = cumsum(g)-g+1
        for(i in seq_along(g)) {
            if(names(g)[i] == "D") {
                o = pos[i]+seq.int(0,length.out=g[i])
                B[o] = "-"
            } else if(names(g)[i] == "I") {
                A = append(A,rep("-", g[i]), pos[i]-1)
                B = append(B,sample(nucs,g[i],prob=p_nucs,replace=TRUE),pos[i]-1)
            }
        }
        am = A != "-"
        bm = B != "-"
        h = ifelse(am, ifelse(bm, "M", "D"), ifelse(bm, "I", "X"))
        h = rle(h)
        stopifnot(h$values == names(g), h$lengths == g)

        sims$mouse = A
        sims$rat = B
        break;
    }
    #modification 
    output = capture.output(write.fasta(sims))
    cat(output,file=paste0(ouDir,basename(input)),append = FALSE,sep='\n')
}

if(!interactive()) {
    ARGS = commandArgs(trailing=TRUE)
    simulate_two_main(input=ARGS[1],ARGS[2])
}
