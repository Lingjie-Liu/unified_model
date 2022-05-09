library(DENR)


### self created example
df1 <- data.frame(seqnames = '1', start = c(4000, 6000),  
                  end= c(10199, 10199),
                  strand = c('+', '+'))
gr1 <- df1 %>% as_granges()
gr1

gr_ds = gr1

#### slide window and tiles
df <- data.frame(seqnames = '1', 
                 start = c(1, 12),  
                 end= c(3, 15),
                 strand = c('+', '+'))
gr <- df %>% as_granges()
plyranges::tile_ranges(gr, width = 1)
plyranges::slide_ranges(gr, width = 1, step = 1)


df <- data.frame(seqnames = '1', 
                 start = c(100, 200, 300, 200, 400, 200, 300),  
                 end= c(100, 200, 300, 200, 400, 200, 300),
                 strand = c('+', '+','-', '+', '+', '+', '-'))

gr <- df %>% plyranges::as_granges() %>% sort()

gr_uniq <- gr %>% unique()
gr_uniq$score <- plyranges::count_overlaps_directed(gr_uniq, gr)






