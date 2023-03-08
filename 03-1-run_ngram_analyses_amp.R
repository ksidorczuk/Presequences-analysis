library(rmarkdown)

render("ngram_analysis_amp.Rmd", output_file = "amp_bigrams_gap0.html", output_dir = "motif-analysis-reports-amp/",
       params = list(k = 2, gaps = 0), fcbf_amp = TRUE)
render("ngram_analysis_amp.Rmd", output_file = "amp_bigrams_gap1.html", output_dir = "motif-analysis-reports-amp/",
       params = list(k = 2, gaps = 1), fcbf_amp = TRUE)
render("ngram_analysis_amp.Rmd", output_file = "amp_bigrams_gap2.html", output_dir = "motif-analysis-reports-amp/",
       params = list(k = 2, gaps = 2), fcbf_amp = TRUE)
render("ngram_analysis_amp.Rmd", output_file = "amp_bigrams_gap3.html", output_dir = "motif-analysis-reports-amp/",
       params = list(k = 2, gaps = 3), fcbf_amp = TRUE)

render("ngram_analysis_amp.Rmd", output_file = "amp_trigrams_gap00.html", output_dir = "motif-analysis-reports-amp/",
       params = list(k = 3, gaps = list(c(0, 0))), fcbf_amp = TRUE)
render("ngram_analysis_amp.Rmd", output_file = "amp_trigrams_gap01.html", output_dir = "motif-analysis-reports-amp/",
       params = list(k = 3, gaps = list(c(0, 1))), fcbf_amp = TRUE)
render("ngram_analysis_amp.Rmd", output_file = "amp_trigrams_gap10.html", output_dir = "motif-analysis-reports-amp/",
       params = list(k = 3, gaps = list(c(1, 0))), fcbf_amp = TRUE)
render("ngram_analysis_amp.Rmd", output_file = "amp_trigrams_gap11.html", output_dir = "motif-analysis-reports-amp/",
       params = list(k = 3, gaps = list(c(1, 1))), fcbf_amp = TRUE)
render("ngram_analysis_amp.Rmd", output_file = "amp_trigrams_gap02.html", output_dir = "motif-analysis-reports-amp/",
       params = list(k = 3, gaps = list(c(0, 2))), fcbf_amp = TRUE)
render("ngram_analysis_amp.Rmd", output_file = "amp_trigrams_gap20.html", output_dir = "motif-analysis-reports-amp/",
       params = list(k = 3, gaps = list(c(2, 0))), fcbf_amp = TRUE)
render("ngram_analysis_amp.Rmd", output_file = "amp_trigrams_gap12.html", output_dir = "motif-analysis-reports-amp/",
       params = list(k = 3, gaps = list(c(1, 2))), fcbf_amp = TRUE)
render("ngram_analysis_amp.Rmd", output_file = "amp_trigrams_gap21.html", output_dir = "motif-analysis-reports-amp/",
       params = list(k = 3, gaps = list(c(2, 1))), fcbf_amp = TRUE)
render("ngram_analysis_amp.Rmd", output_file = "amp_trigrams_gap22.html", output_dir = "motif-analysis-reports-amp/",
       params = list(k = 3, gaps = list(c(2, 2))), fcbf_amp = TRUE)
render("ngram_analysis_amp.Rmd", output_file = "amp_trigrams_gap03.html", output_dir = "motif-analysis-reports-amp/",
       params = list(k = 3, gaps = list(c(0, 3))), fcbf_amp = TRUE)
render("ngram_analysis_amp.Rmd", output_file = "amp_trigrams_gap30.html", output_dir = "motif-analysis-reports-amp/",
       params = list(k = 3, gaps = list(c(3, 0))), fcbf_amp = TRUE)
render("ngram_analysis_amp.Rmd", output_file = "amp_trigrams_gap13.html", output_dir = "motif-analysis-reports-amp/",
       params = list(k = 3, gaps = list(c(1, 3))), fcbf_amp = TRUE)
render("ngram_analysis_amp.Rmd", output_file = "amp_trigrams_gap31.html", output_dir = "motif-analysis-reports-amp/",
       params = list(k = 3, gaps = list(c(3, 1))), fcbf_amp = TRUE)
render("ngram_analysis_amp.Rmd", output_file = "amp_trigrams_gap23.html", output_dir = "motif-analysis-reports-amp/",
       params = list(k = 3, gaps = list(c(2, 3))), fcbf_amp = TRUE)
render("ngram_analysis_amp.Rmd", output_file = "amp_trigrams_gap32.html", output_dir = "motif-analysis-reports-amp/",
       params = list(k = 3, gaps = list(c(3, 2))), fcbf_amp = TRUE)
render("ngram_analysis_amp.Rmd", output_file = "amp_trigrams_gap33.html", output_dir = "motif-analysis-reports-amp/",
       params = list(k = 3, gaps = list(c(3, 3))), fcbf_amp = TRUE)

render("ngram_analysis_amp.Rmd", output_file = "amp_tetragrams_gap000.html", output_dir = "motif-analysis-reports-amp/",
       params = list(k = 4, gaps = list(c(0, 0, 0))), fcbf_amp = TRUE)
render("ngram_analysis_amp.Rmd", output_file = "amp_tetragrams_gap100.html", output_dir = "motif-analysis-reports-amp/",
       params = list(k = 4, gaps = list(c(1, 0, 0))), fcbf_amp = TRUE)
render("ngram_analysis_amp.Rmd", output_file = "amp_tetragrams_gap010.html", output_dir = "motif-analysis-reports-amp/",
       params = list(k = 4, gaps = list(c(0, 1, 0))))
render("ngram_analysis_amp.Rmd", output_file = "amp_tetragrams_gap001.html", output_dir = "motif-analysis-reports-amp/",
       params = list(k = 4, gaps = list(c(0, 0, 1))))
render("ngram_analysis_amp.Rmd", output_file = "amp_tetragrams_gap110.html", output_dir = "motif-analysis-reports-amp/",
       params = list(k = 4, gaps = list(c(1, 1, 0))))
render("ngram_analysis_amp.Rmd", output_file = "amp_tetragrams_gap101.html", output_dir = "motif-analysis-reports-amp/",
       params = list(k = 4, gaps = list(c(1, 0, 1))))
render("ngram_analysis_amp.Rmd", output_file = "amp_tetragrams_gap011.html", output_dir = "motif-analysis-reports-amp/",
       params = list(k = 4, gaps = list(c(0, 1, 1))))
render("ngram_analysis_amp.Rmd", output_file = "amp_tetragrams_gap111.html", output_dir = "motif-analysis-reports-amp/",
       params = list(k = 4, gaps = list(c(1, 1, 1))))
render("ngram_analysis_amp.Rmd", output_file = "amp_tetragrams_gap200.html", output_dir = "motif-analysis-reports-amp/",
       params = list(k = 4, gaps = list(c(2, 0, 0))))
render("ngram_analysis_amp.Rmd", output_file = "amp_tetragrams_gap020.html", output_dir = "motif-analysis-reports-amp/",
       params = list(k = 4, gaps = list(c(0, 2, 0))))
render("ngram_analysis_amp.Rmd", output_file = "amp_tetragrams_gap002.html", output_dir = "motif-analysis-reports-amp/",
       params = list(k = 4, gaps = list(c(0, 0, 2))))
render("ngram_analysis_amp.Rmd", output_file = "amp_tetragrams_gap220.html", output_dir = "motif-analysis-reports-amp/",
       params = list(k = 4, gaps = list(c(2, 2, 0))))
render("ngram_analysis_amp.Rmd", output_file = "amp_tetragrams_gap202.html", output_dir = "motif-analysis-reports-amp/",
       params = list(k = 4, gaps = list(c(2, 0, 2))))
render("ngram_analysis_amp.Rmd", output_file = "amp_tetragrams_gap022.html", output_dir = "motif-analysis-reports-amp/",
       params = list(k = 4, gaps = list(c(0, 2, 2))))
render("ngram_analysis_amp.Rmd", output_file = "amp_tetragrams_gap222.html", output_dir = "motif-analysis-reports-amp/",
       params = list(k = 4, gaps = list(c(2, 2, 2))))
render("ngram_analysis_amp.Rmd", output_file = "amp_tetragrams_gap210.html", output_dir = "motif-analysis-reports-amp/",
       params = list(k = 4, gaps = list(c(2, 1, 0))))
render("ngram_analysis_amp.Rmd", output_file = "amp_tetragrams_gap201.html", output_dir = "motif-analysis-reports-amp/",
       params = list(k = 4, gaps = list(c(2, 0, 1))))
render("ngram_analysis_amp.Rmd", output_file = "amp_tetragrams_gap120.html", output_dir = "motif-analysis-reports-amp/",
       params = list(k = 4, gaps = list(c(1, 2, 0))))
render("ngram_analysis_amp.Rmd", output_file = "amp_tetragrams_gap102.html", output_dir = "motif-analysis-reports-amp/",
       params = list(k = 4, gaps = list(c(1, 0, 2))))
render("ngram_analysis_amp.Rmd", output_file = "amp_tetragrams_gap012.html", output_dir = "motif-analysis-reports-amp/",
       params = list(k = 4, gaps = list(c(0, 1, 2))))
render("ngram_analysis_amp.Rmd", output_file = "amp_tetragrams_gap021.html", output_dir = "motif-analysis-reports-amp/",
       params = list(k = 4, gaps = list(c(0, 2, 1))))
render("ngram_analysis_amp.Rmd", output_file = "amp_tetragrams_gap211.html", output_dir = "motif-analysis-reports-amp/",
       params = list(k = 4, gaps = list(c(2, 1, 1))))
render("ngram_analysis_amp.Rmd", output_file = "amp_tetragrams_gap121.html", output_dir = "motif-analysis-reports-amp/",
       params = list(k = 4, gaps = list(c(1, 2, 1))))
render("ngram_analysis_amp.Rmd", output_file = "amp_tetragrams_gap112.html", output_dir = "motif-analysis-reports-amp/",
       params = list(k = 4, gaps = list(c(1, 1, 2))))
render("ngram_analysis_amp.Rmd", output_file = "amp_tetragrams_gap122.html", output_dir = "motif-analysis-reports-amp/",
       params = list(k = 4, gaps = list(c(1, 2, 2))))
render("ngram_analysis_amp.Rmd", output_file = "amp_tetragrams_gap212.html", output_dir = "motif-analysis-reports-amp/",
       params = list(k = 4, gaps = list(c(2, 1, 2))))
render("ngram_analysis_amp.Rmd", output_file = "amp_tetragrams_gap221.html", output_dir = "motif-analysis-reports-amp/",
       params = list(k = 4, gaps = list(c(2, 2, 1))))

render("ngram_analysis_amp.Rmd", output_file = "amp_pentagrams_gap0000.html", output_dir = "motif-analysis-reports-amp/",
       params = list(k = 5, gaps = list(c(0, 0, 0, 0)), quipt_amp = FALSE))
render("ngram_analysis_amp.Rmd", output_file = "amp_pentagrams_gap1000.html", output_dir = "motif-analysis-reports-amp/",
       params = list(k = 5, gaps = list(c(1, 0, 0, 0)), quipt_amp = FALSE))
render("ngram_analysis_amp.Rmd", output_file = "amp_pentagrams_gap0100.html", output_dir = "motif-analysis-reports-amp/",
       params = list(k = 5, gaps = list(c(0, 1, 0, 0)), quipt_amp = FALSE))
render("ngram_analysis_amp.Rmd", output_file = "amp_pentagrams_gap0010.html", output_dir = "motif-analysis-reports-amp/",
       params = list(k = 5, gaps = list(c(0, 0, 1, 0)), quipt_amp = FALSE))
render("ngram_analysis_amp.Rmd", output_file = "amp_pentagrams_gap0001.html", output_dir = "motif-analysis-reports-amp/",
       params = list(k = 5, gaps = list(c(0, 0, 0, 1)), quipt_amp = FALSE))
render("ngram_analysis_amp.Rmd", output_file = "amp_pentagrams_gap1100.html", output_dir = "motif-analysis-reports-amp/",
       params = list(k = 5, gaps = list(c(1, 1, 0, 0)), quipt_amp = FALSE))
render("ngram_analysis_amp.Rmd", output_file = "amp_pentagrams_gap1010.html", output_dir = "motif-analysis-reports-amp/",
       params = list(k = 5, gaps = list(c(1, 0, 1, 0)), quipt_amp = FALSE))
render("ngram_analysis_amp.Rmd", output_file = "amp_pentagrams_gap1001.html", output_dir = "motif-analysis-reports-amp/",
       params = list(k = 5, gaps = list(c(1, 0, 0, 1)), quipt_amp = FALSE))
render("ngram_analysis_amp.Rmd", output_file = "amp_pentagrams_gap0110.html", output_dir = "motif-analysis-reports-amp/",
       params = list(k = 5, gaps = list(c(0, 1, 1, 0)), quipt_amp = FALSE))
render("ngram_analysis_amp.Rmd", output_file = "amp_pentagrams_gap0101.html", output_dir = "motif-analysis-reports-amp/",
       params = list(k = 5, gaps = list(c(0, 1, 0, 1)), quipt_amp = FALSE))
render("ngram_analysis_amp.Rmd", output_file = "amp_pentagrams_gap0011.html", output_dir = "motif-analysis-reports-amp/",
       params = list(k = 5, gaps = list(c(0, 0, 1, 1)), quipt_amp = FALSE))
render("ngram_analysis_amp.Rmd", output_file = "amp_pentagrams_gap1110.html", output_dir = "motif-analysis-reports-amp/",
       params = list(k = 5, gaps = list(c(1, 1, 1, 0)), quipt_amp = FALSE))
render("ngram_analysis_amp.Rmd", output_file = "amp_pentagrams_gap1011.html", output_dir = "motif-analysis-reports-amp/",
       params = list(k = 5, gaps = list(c(1, 0, 1, 1)), quipt_amp = FALSE))
render("ngram_analysis_amp.Rmd", output_file = "amp_pentagrams_gap1101.html", output_dir = "motif-analysis-reports-amp/",
       params = list(k = 5, gaps = list(c(1, 1, 0, 1)), quipt_amp = FALSE))
render("ngram_analysis_amp.Rmd", output_file = "amp_pentagrams_gap0111.html", output_dir = "motif-analysis-reports-amp/",
       params = list(k = 5, gaps = list(c(0, 1, 1, 1)), quipt_amp = FALSE))
render("ngram_analysis_amp.Rmd", output_file = "amp_pentagrams_gap1111.html", output_dir = "motif-analysis-reports-amp/",
       params = list(k = 5, gaps = list(c(1, 1, 1, 1)), quipt_amp = FALSE))


