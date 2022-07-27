library(rmarkdown)

enc10 <- list("a" = "S", 
              "b" = c("T", "N", "Q"), 
              "c" = "A",
              "d" = c("V", "I", "L"),
              "e" = "M",
              "f" = c("F", "Y", "W"),
              "g" = c("C", "G", "P"),
              "h" = c("D", "E"),
              "i" = "R",
              "j" = c("H", "K"))

enc7 <- list("1" = c("S", "T", "N", "Q"), 
             "3" = c("A", "V", "I", "L"),
             "5" = "M",
             "6" = c("F", "Y", "W"),
             "7" = c("C", "G", "P"),
             "8" = c("D", "E"),
             "9" = c("R", "H", "K"))

enc8 <- list("1" = c("S", "T"),
             "2" = c("N", "Q"), 
             "3" = c("A", "V", "I", "L"),
             "5" = "M",
             "6" = c("F", "Y", "W"),
             "7" = c("C", "G", "P"),
             "8" = c("D", "E"),
             "9" = c("R", "H", "K"))

# encoding enc10
render("ngram_analysis_encoding.Rmd", output_file = "enc10_bigrams_gap0.html", output_path = "motif-analysis-reports-enc/",
       params = list("encoding" = enc10, k = 2, gaps = 0))
render("ngram_analysis_encoding.Rmd", output_file = "enc10_bigrams_gap1.html", output_path = "motif-analysis-reports-enc/",
       params = list("encoding" = enc10, k = 2, gaps = 1))
render("ngram_analysis_encoding.Rmd", output_file = "enc10_bigrams_gap2.html", output_path = "motif-analysis-reports-enc/",
       params = list("encoding" = enc10, k = 2, gaps = 2))
render("ngram_analysis_encoding.Rmd", output_file = "enc10_bigrams_gap3.html", output_path = "motif-analysis-reports-enc/",
       params = list("encoding" = enc10, k = 2, gaps = 3))

render("ngram_analysis_encoding.Rmd", output_file = "enc10_trigrams_gap00.html", output_path = "motif-analysis-reports-enc/",
       params = list("encoding" = enc10, k = 3, gaps = list(c(0, 0))))
render("ngram_analysis_encoding.Rmd", output_file = "enc10_trigrams_gap01.html", output_path = "motif-analysis-reports-enc/",
       params = list("encoding" = enc10, k = 3, gaps = list(c(0, 1))))
render("ngram_analysis_encoding.Rmd", output_file = "enc10_trigrams_gap10.html", output_path = "motif-analysis-reports-enc/",
       params = list("encoding" = enc10, k = 3, gaps = list(c(1, 0))))
render("ngram_analysis_encoding.Rmd", output_file = "enc10_trigrams_gap11.html", output_path = "motif-analysis-reports-enc/",
       params = list("encoding" = enc10, k = 3, gaps = list(c(1, 1))))
render("ngram_analysis_encoding.Rmd", output_file = "enc10_trigrams_gap02.html", output_path = "motif-analysis-reports-enc/",
       params = list("encoding" = enc10, k = 3, gaps = list(c(0, 2))))
render("ngram_analysis_encoding.Rmd", output_file = "enc10_trigrams_gap20.html", output_path = "motif-analysis-reports-enc/",
       params = list("encoding" = enc10, k = 3, gaps = list(c(2, 0))))
render("ngram_analysis_encoding.Rmd", output_file = "enc10_trigrams_gap12.html", output_path = "motif-analysis-reports-enc/",
       params = list("encoding" = enc10, k = 3, gaps = list(c(1, 2))))
render("ngram_analysis_encoding.Rmd", output_file = "enc10_trigrams_gap21.html", output_path = "motif-analysis-reports-enc/",
       params = list("encoding" = enc10, k = 3, gaps = list(c(2, 1))))
render("ngram_analysis_encoding.Rmd", output_file = "enc10_trigrams_gap22.html", output_path = "motif-analysis-reports-enc/",
       params = list("encoding" = enc10, k = 3, gaps = list(c(2, 2))))
render("ngram_analysis_encoding.Rmd", output_file = "enc10_trigrams_gap03.html", output_path = "motif-analysis-reports-enc/",
       params = list("encoding" = enc10, k = 3, gaps = list(c(0, 3))))
render("ngram_analysis_encoding.Rmd", output_file = "enc10_trigrams_gap30.html", output_path = "motif-analysis-reports-enc/",
       params = list("encoding" = enc10, k = 3, gaps = list(c(3, 0))))
render("ngram_analysis_encoding.Rmd", output_file = "enc10_trigrams_gap13.html", output_path = "motif-analysis-reports-enc/",
       params = list("encoding" = enc10, k = 3, gaps = list(c(1, 3))))
render("ngram_analysis_encoding.Rmd", output_file = "enc10_trigrams_gap31.html", output_path = "motif-analysis-reports-enc/",
       params = list("encoding" = enc10, k = 3, gaps = list(c(3, 1))))
render("ngram_analysis_encoding.Rmd", output_file = "enc10_trigrams_gap23.html", output_path = "motif-analysis-reports-enc/",
       params = list("encoding" = enc10, k = 3, gaps = list(c(2, 3))))
render("ngram_analysis_encoding.Rmd", output_file = "enc10_trigrams_gap32.html", output_path = "motif-analysis-reports-enc/",
       params = list("encoding" = enc10, k = 3, gaps = list(c(3, 2))))
render("ngram_analysis_encoding.Rmd", output_file = "enc10_trigrams_gap33.html", output_path = "motif-analysis-reports-enc/",
       params = list("encoding" = enc10, k = 3, gaps = list(c(3, 3))))

render("ngram_analysis_encoding.Rmd", output_file = "enc10_tetragrams_gap000.html", output_dir = "motif-analysis-reports-enc/",
       params = list("encoding" = enc10, k = 4, gaps = list(c(0, 0, 0))))
render("ngram_analysis_encoding.Rmd", output_file = "enc10_tetragrams_gap100.html", output_dir = "motif-analysis-reports-enc/",
       params = list("encoding" = enc10, k = 4, gaps = list(c(1, 0, 0))))
render("ngram_analysis_encoding.Rmd", output_file = "enc10_tetragrams_gap010.html", output_dir = "motif-analysis-reports-enc/",
       params = list("encoding" = enc10, k = 4, gaps = list(c(0, 1, 0))))
render("ngram_analysis_encoding.Rmd", output_file = "enc10_tetragrams_gap001.html", output_dir = "motif-analysis-reports-enc/",
       params = list("encoding" = enc10, k = 4, gaps = list(c(0, 0, 1))))
render("ngram_analysis_encoding.Rmd", output_file = "enc10_tetragrams_gap110.html", output_dir = "motif-analysis-reports-enc/",
       params = list("encoding" = enc10, k = 4, gaps = list(c(1, 1, 0))))
render("ngram_analysis_encoding.Rmd", output_file = "enc10_tetragrams_gap101.html", output_dir = "motif-analysis-reports-enc/",
       params = list("encoding" = enc10, k = 4, gaps = list(c(1, 0, 1))))
render("ngram_analysis_encoding.Rmd", output_file = "enc10_tetragrams_gap011.html", output_dir = "motif-analysis-reports-enc/",
       params = list("encoding" = enc10, k = 4, gaps = list(c(0, 1, 1))))
render("ngram_analysis_encoding.Rmd", output_file = "enc10_tetragrams_gap111.html", output_dir = "motif-analysis-reports-enc/",
       params = list("encoding" = enc10, k = 4, gaps = list(c(1, 1, 1))))
render("ngram_analysis_encoding.Rmd", output_file = "enc10_tetragrams_gap200.html", output_dir = "motif-analysis-reports-enc/",
       params = list("encoding" = enc10, k = 4, gaps = list(c(2, 0, 0))))
render("ngram_analysis_encoding.Rmd", output_file = "enc10_tetragrams_gap020.html", output_dir = "motif-analysis-reports-enc/",
       params = list("encoding" = enc10, k = 4, gaps = list(c(0, 2, 0))))
render("ngram_analysis_encoding.Rmd", output_file = "enc10_tetragrams_gap002.html", output_dir = "motif-analysis-reports-enc/",
       params = list("encoding" = enc10, k = 4, gaps = list(c(0, 0, 2))))
render("ngram_analysis_encoding.Rmd", output_file = "enc10_tetragrams_gap220.html", output_dir = "motif-analysis-reports-enc/",
       params = list("encoding" = enc10, k = 4, gaps = list(c(2, 2, 0))))
render("ngram_analysis_encoding.Rmd", output_file = "enc10_tetragrams_gap202.html", output_dir = "motif-analysis-reports-enc/",
       params = list("encoding" = enc10, k = 4, gaps = list(c(2, 0, 2))))
render("ngram_analysis_encoding.Rmd", output_file = "enc10_tetragrams_gap022.html", output_dir = "motif-analysis-reports-enc/",
       params = list("encoding" = enc10, k = 4, gaps = list(c(0, 2, 2))))
render("ngram_analysis_encoding.Rmd", output_file = "enc10_tetragrams_gap222.html", output_dir = "motif-analysis-reports-enc/",
       params = list("encoding" = enc10, k = 4, gaps = list(c(2, 2, 2))))
render("ngram_analysis_encoding.Rmd", output_file = "enc10_tetragrams_gap210.html", output_dir = "motif-analysis-reports-enc/",
       params = list("encoding" = enc10, k = 4, gaps = list(c(2, 1, 0))))
render("ngram_analysis_encoding.Rmd", output_file = "enc10_tetragrams_gap201.html", output_dir = "motif-analysis-reports-enc/",
       params = list("encoding" = enc10, k = 4, gaps = list(c(2, 0, 1))))
render("ngram_analysis_encoding.Rmd", output_file = "enc10_tetragrams_gap120.html", output_dir = "motif-analysis-reports-enc/",
       params = list("encoding" = enc10, k = 4, gaps = list(c(1, 2, 0))))
render("ngram_analysis_encoding.Rmd", output_file = "enc10_tetragrams_gap102.html", output_dir = "motif-analysis-reports-enc/",
       params = list("encoding" = enc10, k = 4, gaps = list(c(1, 0, 2))))
render("ngram_analysis_encoding.Rmd", output_file = "enc10_tetragrams_gap012.html", output_dir = "motif-analysis-reports-enc/",
       params = list("encoding" = enc10, k = 4, gaps = list(c(0, 1, 2))))
render("ngram_analysis_encoding.Rmd", output_file = "enc10_tetragrams_gap021.html", output_dir = "motif-analysis-reports-enc/",
       params = list("encoding" = enc10, k = 4, gaps = list(c(0, 2, 1))))
render("ngram_analysis_encoding.Rmd", output_file = "enc10_tetragrams_gap211.html", output_dir = "motif-analysis-reports-enc/",
       params = list("encoding" = enc10, k = 4, gaps = list(c(2, 1, 1))))
render("ngram_analysis_encoding.Rmd", output_file = "enc10_tetragrams_gap121.html", output_dir = "motif-analysis-reports-enc/",
       params = list("encoding" = enc10, k = 4, gaps = list(c(1, 2, 1))))
render("ngram_analysis_encoding.Rmd", output_file = "enc10_tetragrams_gap112.html", output_dir = "motif-analysis-reports-enc/",
       params = list("encoding" = enc10, k = 4, gaps = list(c(1, 1, 2))))
render("ngram_analysis_encoding.Rmd", output_file = "enc10_tetragrams_gap122.html", output_dir = "motif-analysis-reports-enc/",
       params = list("encoding" = enc10, k = 4, gaps = list(c(1, 2, 2))))
render("ngram_analysis_encoding.Rmd", output_file = "enc10_tetragrams_gap212.html", output_dir = "motif-analysis-reports-enc/",
       params = list("encoding" = enc10, k = 4, gaps = list(c(2, 1, 2))))
render("ngram_analysis_encoding.Rmd", output_file = "enc10_tetragrams_gap221.html", output_dir = "motif-analysis-reports-enc/",
       params = list("encoding" = enc10, k = 4, gaps = list(c(2, 2, 1))))

render("ngram_analysis_encoding.Rmd", output_file = "enc10_pentagrams_gap0000.html", output_dir = "motif-analysis-reports-enc/",
       params = list(k = 5, gaps = list(c(0, 0, 0, 0))))
render("ngram_analysis_encoding.Rmd", output_file = "enc10_pentagrams_gap1000.html", output_dir = "motif-analysis-reports-enc/",
       params = list(k = 5, gaps = list(c(1, 0, 0, 0))))
render("ngram_analysis_encoding.Rmd", output_file = "enc10_pentagrams_gap0100.html", output_dir = "motif-analysis-reports-enc/",
       params = list(k = 5, gaps = list(c(0, 1, 0, 0))))
render("ngram_analysis_encoding.Rmd", output_file = "enc10_pentagrams_gap0010.html", output_dir = "motif-analysis-reports-enc/",
       params = list(k = 5, gaps = list(c(0, 0, 1, 0))))
render("ngram_analysis_encoding.Rmd", output_file = "enc10_pentagrams_gap0001.html", output_dir = "motif-analysis-reports-enc/",
       params = list(k = 5, gaps = list(c(0, 0, 0, 1))))
render("ngram_analysis_encoding.Rmd", output_file = "enc10_pentagrams_gap1100.html", output_dir = "motif-analysis-reports-enc/",
       params = list(k = 5, gaps = list(c(1, 1, 0, 0))))
render("ngram_analysis_encoding.Rmd", output_file = "enc10_pentagrams_gap1010.html", output_dir = "motif-analysis-reports-enc/",
       params = list(k = 5, gaps = list(c(1, 0, 1, 0))))
render("ngram_analysis_encoding.Rmd", output_file = "enc10_pentagrams_gap1001.html", output_dir = "motif-analysis-reports-enc/",
       params = list(k = 5, gaps = list(c(1, 0, 0, 1))))
render("ngram_analysis_encoding.Rmd", output_file = "enc10_pentagrams_gap0110.html", output_dir = "motif-analysis-reports-enc/",
       params = list(k = 5, gaps = list(c(0, 1, 1, 0))))
render("ngram_analysis_encoding.Rmd", output_file = "enc10_pentagrams_gap0101.html", output_dir = "motif-analysis-reports-enc/",
       params = list(k = 5, gaps = list(c(0, 1, 0, 1))))
render("ngram_analysis_encoding.Rmd", output_file = "enc10_pentagrams_gap0011.html", output_dir = "motif-analysis-reports-enc/",
       params = list(k = 5, gaps = list(c(0, 0, 1, 1))))
render("ngram_analysis_encoding.Rmd", output_file = "enc10_pentagrams_gap1110.html", output_dir = "motif-analysis-reports-enc/",
       params = list(k = 5, gaps = list(c(1, 1, 1, 0))))
render("ngram_analysis_encoding.Rmd", output_file = "enc10_pentagrams_gap1011.html", output_dir = "motif-analysis-reports-enc/",
       params = list(k = 5, gaps = list(c(1, 0, 1, 1))))
render("ngram_analysis_encoding.Rmd", output_file = "enc10_pentagrams_gap1101.html", output_dir = "motif-analysis-reports-enc/",
       params = list(k = 5, gaps = list(c(1, 1, 0, 1))))
render("ngram_analysis_encoding.Rmd", output_file = "enc10_pentagrams_gap0111.html", output_dir = "motif-analysis-reports-enc/",
       params = list(k = 5, gaps = list(c(0, 1, 1, 1))))
render("ngram_analysis_encoding.Rmd", output_file = "enc10_pentagrams_gap1111.html", output_dir = "motif-analysis-reports-enc/",
       params = list(k = 5, gaps = list(c(1, 1, 1, 1))))


# encoding enc8
render("ngram_analysis_encoding.Rmd", output_file = "enc8_bigrams_gap0.html", output_path = "motif-analysis-reports-enc/",
       params = list("encoding" = enc8, k = 2, gaps = 0))
render("ngram_analysis_encoding.Rmd", output_file = "enc8_bigrams_gap1.html", output_path = "motif-analysis-reports-enc/",
       params = list("encoding" = enc8, k = 2, gaps = 1))
render("ngram_analysis_encoding.Rmd", output_file = "enc8_bigrams_gap2.html", output_path = "motif-analysis-reports-enc/",
       params = list("encoding" = enc8, k = 2, gaps = 2))
render("ngram_analysis_encoding.Rmd", output_file = "enc8_bigrams_gap3.html", output_path = "motif-analysis-reports-enc/",
       params = list("encoding" = enc8, k = 2, gaps = 3))

render("ngram_analysis_encoding.Rmd", output_file = "enc8_trigrams_gap00.html", output_path = "motif-analysis-reports-enc/",
       params = list("encoding" = enc8, k = 3, gaps = list(c(0, 0))))
render("ngram_analysis_encoding.Rmd", output_file = "enc8_trigrams_gap01.html", output_path = "motif-analysis-reports-enc/",
       params = list("encoding" = enc8, k = 3, gaps = list(c(0, 1))))
render("ngram_analysis_encoding.Rmd", output_file = "enc8_trigrams_gap10.html", output_path = "motif-analysis-reports-enc/",
       params = list("encoding" = enc8, k = 3, gaps = list(c(1, 0))))
render("ngram_analysis_encoding.Rmd", output_file = "enc8_trigrams_gap11.html", output_path = "motif-analysis-reports-enc/",
       params = list("encoding" = enc8, k = 3, gaps = list(c(1, 1))))
render("ngram_analysis_encoding.Rmd", output_file = "enc8_trigrams_gap02.html", output_path = "motif-analysis-reports-enc/",
       params = list("encoding" = enc8, k = 3, gaps = list(c(0, 2))))
render("ngram_analysis_encoding.Rmd", output_file = "enc8_trigrams_gap20.html", output_path = "motif-analysis-reports-enc/",
       params = list("encoding" = enc8, k = 3, gaps = list(c(2, 0))))
render("ngram_analysis_encoding.Rmd", output_file = "enc8_trigrams_gap12.html", output_path = "motif-analysis-reports-enc/",
       params = list("encoding" = enc8, k = 3, gaps = list(c(1, 2))))
render("ngram_analysis_encoding.Rmd", output_file = "enc8_trigrams_gap21.html", output_path = "motif-analysis-reports-enc/",
       params = list("encoding" = enc8, k = 3, gaps = list(c(2, 1))))
render("ngram_analysis_encoding.Rmd", output_file = "enc8_trigrams_gap22.html", output_path = "motif-analysis-reports-enc/",
       params = list("encoding" = enc8, k = 3, gaps = list(c(2, 2))))
render("ngram_analysis_encoding.Rmd", output_file = "enc8_trigrams_gap03.html", output_path = "motif-analysis-reports-enc/",
       params = list("encoding" = enc8, k = 3, gaps = list(c(0, 3))))
render("ngram_analysis_encoding.Rmd", output_file = "enc8_trigrams_gap30.html", output_path = "motif-analysis-reports-enc/",
       params = list("encoding" = enc8, k = 3, gaps = list(c(3, 0))))
render("ngram_analysis_encoding.Rmd", output_file = "enc8_trigrams_gap13.html", output_path = "motif-analysis-reports-enc/",
       params = list("encoding" = enc8, k = 3, gaps = list(c(1, 3))))
render("ngram_analysis_encoding.Rmd", output_file = "enc8_trigrams_gap31.html", output_path = "motif-analysis-reports-enc/",
       params = list("encoding" = enc8, k = 3, gaps = list(c(3, 1))))
render("ngram_analysis_encoding.Rmd", output_file = "enc8_trigrams_gap23.html", output_path = "motif-analysis-reports-enc/",
       params = list("encoding" = enc8, k = 3, gaps = list(c(2, 3))))
render("ngram_analysis_encoding.Rmd", output_file = "enc8_trigrams_gap32.html", output_path = "motif-analysis-reports-enc/",
       params = list("encoding" = enc8, k = 3, gaps = list(c(3, 2))))
render("ngram_analysis_encoding.Rmd", output_file = "enc8_trigrams_gap33.html", output_path = "motif-analysis-reports-enc/",
       params = list("encoding" = enc8, k = 3, gaps = list(c(3, 3))))

render("ngram_analysis_encoding.Rmd", output_file = "enc8_tetragrams_gap000.html", output_dir = "motif-analysis-reports-enc/",
       params = list("encoding" = enc8, k = 4, gaps = list(c(0, 0, 0))))
render("ngram_analysis_encoding.Rmd", output_file = "enc8_tetragrams_gap100.html", output_dir = "motif-analysis-reports-enc/",
       params = list("encoding" = enc8, k = 4, gaps = list(c(1, 0, 0))))
render("ngram_analysis_encoding.Rmd", output_file = "enc8_tetragrams_gap010.html", output_dir = "motif-analysis-reports-enc/",
       params = list("encoding" = enc8, k = 4, gaps = list(c(0, 1, 0))))
render("ngram_analysis_encoding.Rmd", output_file = "enc8_tetragrams_gap001.html", output_dir = "motif-analysis-reports-enc/",
       params = list("encoding" = enc8, k = 4, gaps = list(c(0, 0, 1))))
render("ngram_analysis_encoding.Rmd", output_file = "enc8_tetragrams_gap110.html", output_dir = "motif-analysis-reports-enc/",
       params = list("encoding" = enc8, k = 4, gaps = list(c(1, 1, 0))))
render("ngram_analysis_encoding.Rmd", output_file = "enc8_tetragrams_gap101.html", output_dir = "motif-analysis-reports-enc/",
       params = list("encoding" = enc8, k = 4, gaps = list(c(1, 0, 1))))
render("ngram_analysis_encoding.Rmd", output_file = "enc8_tetragrams_gap011.html", output_dir = "motif-analysis-reports-enc/",
       params = list("encoding" = enc8, k = 4, gaps = list(c(0, 1, 1))))
render("ngram_analysis_encoding.Rmd", output_file = "enc8_tetragrams_gap111.html", output_dir = "motif-analysis-reports-enc/",
       params = list("encoding" = enc8, k = 4, gaps = list(c(1, 1, 1))))
render("ngram_analysis_encoding.Rmd", output_file = "enc8_tetragrams_gap200.html", output_dir = "motif-analysis-reports-enc/",
       params = list("encoding" = enc8, k = 4, gaps = list(c(2, 0, 0))))
render("ngram_analysis_encoding.Rmd", output_file = "enc8_tetragrams_gap020.html", output_dir = "motif-analysis-reports-enc/",
       params = list("encoding" = enc8, k = 4, gaps = list(c(0, 2, 0))))
render("ngram_analysis_encoding.Rmd", output_file = "enc8_tetragrams_gap002.html", output_dir = "motif-analysis-reports-enc/",
       params = list("encoding" = enc8, k = 4, gaps = list(c(0, 0, 2))))
render("ngram_analysis_encoding.Rmd", output_file = "enc8_tetragrams_gap220.html", output_dir = "motif-analysis-reports-enc/",
       params = list("encoding" = enc8, k = 4, gaps = list(c(2, 2, 0))))
render("ngram_analysis_encoding.Rmd", output_file = "enc8_tetragrams_gap202.html", output_dir = "motif-analysis-reports-enc/",
       params = list("encoding" = enc8, k = 4, gaps = list(c(2, 0, 2))))
render("ngram_analysis_encoding.Rmd", output_file = "enc8_tetragrams_gap022.html", output_dir = "motif-analysis-reports-enc/",
       params = list("encoding" = enc8, k = 4, gaps = list(c(0, 2, 2))))
render("ngram_analysis_encoding.Rmd", output_file = "enc8_tetragrams_gap222.html", output_dir = "motif-analysis-reports-enc/",
       params = list("encoding" = enc8, k = 4, gaps = list(c(2, 2, 2))))
render("ngram_analysis_encoding.Rmd", output_file = "enc8_tetragrams_gap210.html", output_dir = "motif-analysis-reports-enc/",
       params = list("encoding" = enc8, k = 4, gaps = list(c(2, 1, 0))))
render("ngram_analysis_encoding.Rmd", output_file = "enc8_tetragrams_gap201.html", output_dir = "motif-analysis-reports-enc/",
       params = list("encoding" = enc8, k = 4, gaps = list(c(2, 0, 1))))
render("ngram_analysis_encoding.Rmd", output_file = "enc8_tetragrams_gap120.html", output_dir = "motif-analysis-reports-enc/",
       params = list("encoding" = enc8, k = 4, gaps = list(c(1, 2, 0))))
render("ngram_analysis_encoding.Rmd", output_file = "enc8_tetragrams_gap102.html", output_dir = "motif-analysis-reports-enc/",
       params = list("encoding" = enc8, k = 4, gaps = list(c(1, 0, 2))))
render("ngram_analysis_encoding.Rmd", output_file = "enc8_tetragrams_gap012.html", output_dir = "motif-analysis-reports-enc/",
       params = list("encoding" = enc8, k = 4, gaps = list(c(0, 1, 2))))
render("ngram_analysis_encoding.Rmd", output_file = "enc8_tetragrams_gap021.html", output_dir = "motif-analysis-reports-enc/",
       params = list("encoding" = enc8, k = 4, gaps = list(c(0, 2, 1))))
render("ngram_analysis_encoding.Rmd", output_file = "enc8_tetragrams_gap211.html", output_dir = "motif-analysis-reports-enc/",
       params = list("encoding" = enc8, k = 4, gaps = list(c(2, 1, 1))))
render("ngram_analysis_encoding.Rmd", output_file = "enc8_tetragrams_gap121.html", output_dir = "motif-analysis-reports-enc/",
       params = list("encoding" = enc8, k = 4, gaps = list(c(1, 2, 1))))
render("ngram_analysis_encoding.Rmd", output_file = "enc8_tetragrams_gap112.html", output_dir = "motif-analysis-reports-enc/",
       params = list("encoding" = enc8, k = 4, gaps = list(c(1, 1, 2))))
render("ngram_analysis_encoding.Rmd", output_file = "enc8_tetragrams_gap122.html", output_dir = "motif-analysis-reports-enc/",
       params = list("encoding" = enc8, k = 4, gaps = list(c(1, 2, 2))))
render("ngram_analysis_encoding.Rmd", output_file = "enc8_tetragrams_gap212.html", output_dir = "motif-analysis-reports-enc/",
       params = list("encoding" = enc8, k = 4, gaps = list(c(2, 1, 2))))
render("ngram_analysis_encoding.Rmd", output_file = "enc8_tetragrams_gap221.html", output_dir = "motif-analysis-reports-enc/",
       params = list("encoding" = enc8, k = 4, gaps = list(c(2, 2, 1))))

render("ngram_analysis_encoding.Rmd", output_file = "enc8_pentagrams_gap0000.html", output_dir = "motif-analysis-reports-enc/",
       params = list(k = 5, gaps = list(c(0, 0, 0, 0))))
render("ngram_analysis_encoding.Rmd", output_file = "enc8_pentagrams_gap1000.html", output_dir = "motif-analysis-reports-enc/",
       params = list(k = 5, gaps = list(c(1, 0, 0, 0))))
render("ngram_analysis_encoding.Rmd", output_file = "enc8_pentagrams_gap0100.html", output_dir = "motif-analysis-reports-enc/",
       params = list(k = 5, gaps = list(c(0, 1, 0, 0))))
render("ngram_analysis_encoding.Rmd", output_file = "enc8_pentagrams_gap0010.html", output_dir = "motif-analysis-reports-enc/",
       params = list(k = 5, gaps = list(c(0, 0, 1, 0))))
render("ngram_analysis_encoding.Rmd", output_file = "enc8_pentagrams_gap0001.html", output_dir = "motif-analysis-reports-enc/",
       params = list(k = 5, gaps = list(c(0, 0, 0, 1))))
render("ngram_analysis_encoding.Rmd", output_file = "enc8_pentagrams_gap1100.html", output_dir = "motif-analysis-reports-enc/",
       params = list(k = 5, gaps = list(c(1, 1, 0, 0))))
render("ngram_analysis_encoding.Rmd", output_file = "enc8_pentagrams_gap1010.html", output_dir = "motif-analysis-reports-enc/",
       params = list(k = 5, gaps = list(c(1, 0, 1, 0))))
render("ngram_analysis_encoding.Rmd", output_file = "enc8_pentagrams_gap1001.html", output_dir = "motif-analysis-reports-enc/",
       params = list(k = 5, gaps = list(c(1, 0, 0, 1))))
render("ngram_analysis_encoding.Rmd", output_file = "enc8_pentagrams_gap0110.html", output_dir = "motif-analysis-reports-enc/",
       params = list(k = 5, gaps = list(c(0, 1, 1, 0))))
render("ngram_analysis_encoding.Rmd", output_file = "enc8_pentagrams_gap0101.html", output_dir = "motif-analysis-reports-enc/",
       params = list(k = 5, gaps = list(c(0, 1, 0, 1))))
render("ngram_analysis_encoding.Rmd", output_file = "enc8_pentagrams_gap0011.html", output_dir = "motif-analysis-reports-enc/",
       params = list(k = 5, gaps = list(c(0, 0, 1, 1))))
render("ngram_analysis_encoding.Rmd", output_file = "enc8_pentagrams_gap1110.html", output_dir = "motif-analysis-reports-enc/",
       params = list(k = 5, gaps = list(c(1, 1, 1, 0))))
render("ngram_analysis_encoding.Rmd", output_file = "enc8_pentagrams_gap1011.html", output_dir = "motif-analysis-reports-enc/",
       params = list(k = 5, gaps = list(c(1, 0, 1, 1))))
render("ngram_analysis_encoding.Rmd", output_file = "enc8_pentagrams_gap1101.html", output_dir = "motif-analysis-reports-enc/",
       params = list(k = 5, gaps = list(c(1, 1, 0, 1))))
render("ngram_analysis_encoding.Rmd", output_file = "enc8_pentagrams_gap0111.html", output_dir = "motif-analysis-reports-enc/",
       params = list(k = 5, gaps = list(c(0, 1, 1, 1))))
render("ngram_analysis_encoding.Rmd", output_file = "enc8_pentagrams_gap1111.html", output_dir = "motif-analysis-reports-enc/",
       params = list(k = 5, gaps = list(c(1, 1, 1, 1))))


# encoding enc7
render("ngram_analysis_encoding.Rmd", output_file = "enc7_bigrams_gap0.html", output_path = "motif-analysis-reports-enc/",
       params = list("encoding" = enc7, k = 2, gaps = 0))
render("ngram_analysis_encoding.Rmd", output_file = "enc7_bigrams_gap1.html", output_path = "motif-analysis-reports-enc/",
       params = list("encoding" = enc7, k = 2, gaps = 1))
render("ngram_analysis_encoding.Rmd", output_file = "enc7_bigrams_gap2.html", output_path = "motif-analysis-reports-enc/",
       params = list("encoding" = enc7, k = 2, gaps = 2))
render("ngram_analysis_encoding.Rmd", output_file = "enc7_bigrams_gap3.html", output_path = "motif-analysis-reports-enc/",
       params = list("encoding" = enc7, k = 2, gaps = 3))

render("ngram_analysis_encoding.Rmd", output_file = "enc7_trigrams_gap00.html", output_path = "motif-analysis-reports-enc/",
       params = list("encoding" = enc7, k = 3, gaps = list(c(0, 0))))
render("ngram_analysis_encoding.Rmd", output_file = "enc7_trigrams_gap01.html", output_path = "motif-analysis-reports-enc/",
       params = list("encoding" = enc7, k = 3, gaps = list(c(0, 1))))
render("ngram_analysis_encoding.Rmd", output_file = "enc7_trigrams_gap10.html", output_path = "motif-analysis-reports-enc/",
       params = list("encoding" = enc7, k = 3, gaps = list(c(1, 0))))
render("ngram_analysis_encoding.Rmd", output_file = "enc7_trigrams_gap11.html", output_path = "motif-analysis-reports-enc/",
       params = list("encoding" = enc7, k = 3, gaps = list(c(1, 1))))
render("ngram_analysis_encoding.Rmd", output_file = "enc7_trigrams_gap02.html", output_path = "motif-analysis-reports-enc/",
       params = list("encoding" = enc7, k = 3, gaps = list(c(0, 2))))
render("ngram_analysis_encoding.Rmd", output_file = "enc7_trigrams_gap20.html", output_path = "motif-analysis-reports-enc/",
       params = list("encoding" = enc7, k = 3, gaps = list(c(2, 0))))
render("ngram_analysis_encoding.Rmd", output_file = "enc7_trigrams_gap12.html", output_path = "motif-analysis-reports-enc/",
       params = list("encoding" = enc7, k = 3, gaps = list(c(1, 2))))
render("ngram_analysis_encoding.Rmd", output_file = "enc7_trigrams_gap21.html", output_path = "motif-analysis-reports-enc/",
       params = list("encoding" = enc7, k = 3, gaps = list(c(2, 1))))
render("ngram_analysis_encoding.Rmd", output_file = "enc7_trigrams_gap22.html", output_path = "motif-analysis-reports-enc/",
       params = list("encoding" = enc7, k = 3, gaps = list(c(2, 2))))
render("ngram_analysis_encoding.Rmd", output_file = "enc7_trigrams_gap03.html", output_path = "motif-analysis-reports-enc/",
       params = list("encoding" = enc7, k = 3, gaps = list(c(0, 3))))
render("ngram_analysis_encoding.Rmd", output_file = "enc7_trigrams_gap30.html", output_path = "motif-analysis-reports-enc/",
       params = list("encoding" = enc7, k = 3, gaps = list(c(3, 0))))
render("ngram_analysis_encoding.Rmd", output_file = "enc7_trigrams_gap13.html", output_path = "motif-analysis-reports-enc/",
       params = list("encoding" = enc7, k = 3, gaps = list(c(1, 3))))
render("ngram_analysis_encoding.Rmd", output_file = "enc7_trigrams_gap31.html", output_path = "motif-analysis-reports-enc/",
       params = list("encoding" = enc7, k = 3, gaps = list(c(3, 1))))
render("ngram_analysis_encoding.Rmd", output_file = "enc7_trigrams_gap23.html", output_path = "motif-analysis-reports-enc/",
       params = list("encoding" = enc7, k = 3, gaps = list(c(2, 3))))
render("ngram_analysis_encoding.Rmd", output_file = "enc7_trigrams_gap32.html", output_path = "motif-analysis-reports-enc/",
       params = list("encoding" = enc7, k = 3, gaps = list(c(3, 2))))
render("ngram_analysis_encoding.Rmd", output_file = "enc7_trigrams_gap33.html", output_path = "motif-analysis-reports-enc/",
       params = list("encoding" = enc7, k = 3, gaps = list(c(3, 3))))

render("ngram_analysis_encoding.Rmd", output_file = "enc7_tetragrams_gap000.html", output_dir = "motif-analysis-reports-enc/",
       params = list("encoding" = enc7, k = 4, gaps = list(c(0, 0, 0))))
render("ngram_analysis_encoding.Rmd", output_file = "enc7_tetragrams_gap100.html", output_dir = "motif-analysis-reports-enc/",
       params = list("encoding" = enc7, k = 4, gaps = list(c(1, 0, 0))))
render("ngram_analysis_encoding.Rmd", output_file = "enc7_tetragrams_gap010.html", output_dir = "motif-analysis-reports-enc/",
       params = list("encoding" = enc7, k = 4, gaps = list(c(0, 1, 0))))
render("ngram_analysis_encoding.Rmd", output_file = "enc7_tetragrams_gap001.html", output_dir = "motif-analysis-reports-enc/",
       params = list("encoding" = enc7, k = 4, gaps = list(c(0, 0, 1))))
render("ngram_analysis_encoding.Rmd", output_file = "enc7_tetragrams_gap110.html", output_dir = "motif-analysis-reports-enc/",
       params = list("encoding" = enc7, k = 4, gaps = list(c(1, 1, 0))))
render("ngram_analysis_encoding.Rmd", output_file = "enc7_tetragrams_gap101.html", output_dir = "motif-analysis-reports-enc/",
       params = list("encoding" = enc7, k = 4, gaps = list(c(1, 0, 1))))
render("ngram_analysis_encoding.Rmd", output_file = "enc7_tetragrams_gap011.html", output_dir = "motif-analysis-reports-enc/",
       params = list("encoding" = enc7, k = 4, gaps = list(c(0, 1, 1))))
render("ngram_analysis_encoding.Rmd", output_file = "enc7_tetragrams_gap111.html", output_dir = "motif-analysis-reports-enc/",
       params = list("encoding" = enc7, k = 4, gaps = list(c(1, 1, 1))))
render("ngram_analysis_encoding.Rmd", output_file = "enc7_tetragrams_gap200.html", output_dir = "motif-analysis-reports-enc/",
       params = list("encoding" = enc7, k = 4, gaps = list(c(2, 0, 0))))
render("ngram_analysis_encoding.Rmd", output_file = "enc7_tetragrams_gap020.html", output_dir = "motif-analysis-reports-enc/",
       params = list("encoding" = enc7, k = 4, gaps = list(c(0, 2, 0))))
render("ngram_analysis_encoding.Rmd", output_file = "enc7_tetragrams_gap002.html", output_dir = "motif-analysis-reports-enc/",
       params = list("encoding" = enc7, k = 4, gaps = list(c(0, 0, 2))))
render("ngram_analysis_encoding.Rmd", output_file = "enc7_tetragrams_gap220.html", output_dir = "motif-analysis-reports-enc/",
       params = list("encoding" = enc7, k = 4, gaps = list(c(2, 2, 0))))
render("ngram_analysis_encoding.Rmd", output_file = "enc7_tetragrams_gap202.html", output_dir = "motif-analysis-reports-enc/",
       params = list("encoding" = enc7, k = 4, gaps = list(c(2, 0, 2))))
render("ngram_analysis_encoding.Rmd", output_file = "enc7_tetragrams_gap022.html", output_dir = "motif-analysis-reports-enc/",
       params = list("encoding" = enc7, k = 4, gaps = list(c(0, 2, 2))))
render("ngram_analysis_encoding.Rmd", output_file = "enc7_tetragrams_gap222.html", output_dir = "motif-analysis-reports-enc/",
       params = list("encoding" = enc7, k = 4, gaps = list(c(2, 2, 2))))
render("ngram_analysis_encoding.Rmd", output_file = "enc7_tetragrams_gap210.html", output_dir = "motif-analysis-reports-enc/",
       params = list("encoding" = enc7, k = 4, gaps = list(c(2, 1, 0))))
render("ngram_analysis_encoding.Rmd", output_file = "enc7_tetragrams_gap201.html", output_dir = "motif-analysis-reports-enc/",
       params = list("encoding" = enc7, k = 4, gaps = list(c(2, 0, 1))))
render("ngram_analysis_encoding.Rmd", output_file = "enc7_tetragrams_gap120.html", output_dir = "motif-analysis-reports-enc/",
       params = list("encoding" = enc7, k = 4, gaps = list(c(1, 2, 0))))
render("ngram_analysis_encoding.Rmd", output_file = "enc7_tetragrams_gap102.html", output_dir = "motif-analysis-reports-enc/",
       params = list("encoding" = enc7, k = 4, gaps = list(c(1, 0, 2))))
render("ngram_analysis_encoding.Rmd", output_file = "enc7_tetragrams_gap012.html", output_dir = "motif-analysis-reports-enc/",
       params = list("encoding" = enc7, k = 4, gaps = list(c(0, 1, 2))))
render("ngram_analysis_encoding.Rmd", output_file = "enc7_tetragrams_gap021.html", output_dir = "motif-analysis-reports-enc/",
       params = list("encoding" = enc7, k = 4, gaps = list(c(0, 2, 1))))
render("ngram_analysis_encoding.Rmd", output_file = "enc7_tetragrams_gap211.html", output_dir = "motif-analysis-reports-enc/",
       params = list("encoding" = enc7, k = 4, gaps = list(c(2, 1, 1))))
render("ngram_analysis_encoding.Rmd", output_file = "enc7_tetragrams_gap121.html", output_dir = "motif-analysis-reports-enc/",
       params = list("encoding" = enc7, k = 4, gaps = list(c(1, 2, 1))))
render("ngram_analysis_encoding.Rmd", output_file = "enc7_tetragrams_gap112.html", output_dir = "motif-analysis-reports-enc/",
       params = list("encoding" = enc7, k = 4, gaps = list(c(1, 1, 2))))
render("ngram_analysis_encoding.Rmd", output_file = "enc7_tetragrams_gap122.html", output_dir = "motif-analysis-reports-enc/",
       params = list("encoding" = enc7, k = 4, gaps = list(c(1, 2, 2))))
render("ngram_analysis_encoding.Rmd", output_file = "enc7_tetragrams_gap212.html", output_dir = "motif-analysis-reports-enc/",
       params = list("encoding" = enc7, k = 4, gaps = list(c(2, 1, 2))))
render("ngram_analysis_encoding.Rmd", output_file = "enc7_tetragrams_gap221.html", output_dir = "motif-analysis-reports-enc/",
       params = list("encoding" = enc7, k = 4, gaps = list(c(2, 2, 1))))

render("ngram_analysis_encoding.Rmd", output_file = "enc7_pentagrams_gap0000.html", output_dir = "motif-analysis-reports-enc/",
       params = list(k = 5, gaps = list(c(0, 0, 0, 0))))
render("ngram_analysis_encoding.Rmd", output_file = "enc7_pentagrams_gap1000.html", output_dir = "motif-analysis-reports-enc/",
       params = list(k = 5, gaps = list(c(1, 0, 0, 0))))
render("ngram_analysis_encoding.Rmd", output_file = "enc7_pentagrams_gap0100.html", output_dir = "motif-analysis-reports-enc/",
       params = list(k = 5, gaps = list(c(0, 1, 0, 0))))
render("ngram_analysis_encoding.Rmd", output_file = "enc7_pentagrams_gap0010.html", output_dir = "motif-analysis-reports-enc/",
       params = list(k = 5, gaps = list(c(0, 0, 1, 0))))
render("ngram_analysis_encoding.Rmd", output_file = "enc7_pentagrams_gap0001.html", output_dir = "motif-analysis-reports-enc/",
       params = list(k = 5, gaps = list(c(0, 0, 0, 1))))
render("ngram_analysis_encoding.Rmd", output_file = "enc7_pentagrams_gap1100.html", output_dir = "motif-analysis-reports-enc/",
       params = list(k = 5, gaps = list(c(1, 1, 0, 0))))
render("ngram_analysis_encoding.Rmd", output_file = "enc7_pentagrams_gap1010.html", output_dir = "motif-analysis-reports-enc/",
       params = list(k = 5, gaps = list(c(1, 0, 1, 0))))
render("ngram_analysis_encoding.Rmd", output_file = "enc7_pentagrams_gap1001.html", output_dir = "motif-analysis-reports-enc/",
       params = list(k = 5, gaps = list(c(1, 0, 0, 1))))
render("ngram_analysis_encoding.Rmd", output_file = "enc7_pentagrams_gap0110.html", output_dir = "motif-analysis-reports-enc/",
       params = list(k = 5, gaps = list(c(0, 1, 1, 0))))
render("ngram_analysis_encoding.Rmd", output_file = "enc7_pentagrams_gap0101.html", output_dir = "motif-analysis-reports-enc/",
       params = list(k = 5, gaps = list(c(0, 1, 0, 1))))
render("ngram_analysis_encoding.Rmd", output_file = "enc7_pentagrams_gap0011.html", output_dir = "motif-analysis-reports-enc/",
       params = list(k = 5, gaps = list(c(0, 0, 1, 1))))
render("ngram_analysis_encoding.Rmd", output_file = "enc7_pentagrams_gap1110.html", output_dir = "motif-analysis-reports-enc/",
       params = list(k = 5, gaps = list(c(1, 1, 1, 0))))
render("ngram_analysis_encoding.Rmd", output_file = "enc7_pentagrams_gap1011.html", output_dir = "motif-analysis-reports-enc/",
       params = list(k = 5, gaps = list(c(1, 0, 1, 1))))
render("ngram_analysis_encoding.Rmd", output_file = "enc7_pentagrams_gap1101.html", output_dir = "motif-analysis-reports-enc/",
       params = list(k = 5, gaps = list(c(1, 1, 0, 1))))
render("ngram_analysis_encoding.Rmd", output_file = "enc7_pentagrams_gap0111.html", output_dir = "motif-analysis-reports-enc/",
       params = list(k = 5, gaps = list(c(0, 1, 1, 1))))
render("ngram_analysis_encoding.Rmd", output_file = "enc7_pentagrams_gap1111.html", output_dir = "motif-analysis-reports-enc/",
       params = list(k = 5, gaps = list(c(1, 1, 1, 1))))