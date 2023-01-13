library(rmarkdown)

render("ngram_analysis_counts.Rmd", output_file = "counts_bigrams_without_gaps.html", output_dir = "motif-analysis-reports-counts/",
       params = list(k = 2, gaps = 0))
render("ngram_analysis_counts.Rmd", output_file = "counts_bigrams_with_gaps.html", output_dir = "motif-analysis-reports-counts/",
       params = list(k = c(2, 2, 2), gaps = c(0, 1, 2)))

render("ngram_analysis_counts.Rmd", output_file = "counts_trigrams_without_gaps.html", output_dir = "motif-analysis-reports-counts/",
       params = list(k = 3, gaps = list(c(0, 0))))
render("ngram_analysis_counts.Rmd", output_file = "counts_trigrams_with_gaps.html", output_dir = "motif-analysis-reports-counts/",
       params = list(k = rep(c(3), 16), gaps = list(c(0, 0), c(0, 1), c(1, 0), c(1, 1), 
                                                    c(0, 2), c(2, 0), c(1, 2), c(2, 1), 
                                                    c(2, 2), c(0, 3), c(3, 0), c(1, 3), 
                                                    c(3, 1), c(2, 3), c(3, 2), c(3, 3))))

render("ngram_analysis_counts.Rmd", output_file = "counts_tetragrams_without_gaps.html", output_dir = "motif-analysis-reports-counts/",
       params = list(k = 4, gaps = list(c(0, 0, 0))))
render("ngram_analysis_counts.Rmd", output_file = "counts_tetragrams_with_gaps.html", output_dir = "motif-analysis-reports-counts/",
       params = list(k = rep(c(4), 27), gaps = list(c(0, 0, 0), c(1, 0, 0), c(0, 1, 0), c(0, 0, 1), 
                                                    c(1, 1, 0), c(1, 0, 1), c(0, 1, 1), c(1, 1, 1), 
                                                    c(2, 0, 0), c(0, 2, 0), c(0, 0, 2), c(2, 2, 0), 
                                                    c(2, 0, 2), c(0, 2, 2), c(2, 2, 2), c(2, 1, 0), 
                                                    c(2, 0, 1), c(1, 2, 0), c(1, 0, 2), c(0, 1, 2), 
                                                    c(0, 2, 1), c(2, 1, 1), c(1, 2, 1), c(1, 1, 2), 
                                                    c(1, 2, 2), c(2, 1, 2), c(2, 2, 1))))

render("ngram_analysis_counts.Rmd", output_file = "counts_pentagrams_without_gaps.html", output_dir = "motif-analysis-reports-counts/",
       params = list(k = 5, gaps = list(c(0, 0, 0, 0))))
render("ngram_analysis_counts.Rmd", output_file = "counts_pentagrams_with_gaps.html", output_dir = "motif-analysis-reports-counts/",
       params = list(k = rep(c(5), 16), gaps = list(c(0, 0, 0, 0), c(1, 0, 0, 0), c(0, 1, 0, 0), c(0, 0, 1, 0), 
                                                    c(0, 0, 0, 1), c(1, 1, 0, 0), c(1, 0, 1, 0), c(1, 0, 0, 1), 
                                                    c(0, 1, 1, 0), c(0, 1, 0, 1), c(0, 0, 1, 1), c(1, 1, 1, 0), 
                                                    c(1, 0, 1, 1), c(1, 1, 0, 1), c(0, 1, 1, 1), c(1, 1, 1, 1))))
