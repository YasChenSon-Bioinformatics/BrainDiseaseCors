
library(network)
library(sna)
n <- network(rgraph(10, tprob = 0.2), directed = FALSE)
n %v% "family" <- sample(letters[1:3], 10, replace = TRUE)
n %v% "family"
n %v% "importance" <- sample(1:3, 10, replace = TRUE)
e <- network.edgecount(n)
set.edge.attribute(n, "type", sample(letters[24:26], e, replace = TRUE))
set.edge.attribute(n, "day", sample(1:3, e, replace = TRUE))
ggnetwork(n, layout = "fruchtermanreingold", cell.jitter = 0.75)
ggnetwork(n, layout = "target", niter = 100)

rownames(deg_matrix)

ggnetwork(dnet, layout = fixed_deg_layout)

    ggplot(dnet, aes(x = x, y = y, xend = xend, yend = yend)) +
    geom_edges(color = "grey50", aes(size=n_DEG_up+n_DEG_lo), show.legend = FALSE ) +
    geom_edgetext_repel(aes(label = n_DEG_up), color = "black", size = 7,
                        fill = "#E6AAAA", box.padding = unit(1.5, "lines") ) +
    geom_edgetext_repel(aes(label = n_DEG_lo), color = "black", size = 7,
                        fill = "skyblue", box.padding = unit(.5, "lines")) +
    geom_nodes(color = "white", size = 17) +
    geom_nodetext(aes(label = paste0("GDS\n",substr(vertex.names,4,100)) ),
                  fontface = "bold", color="black", size=8, lineheight=.8) +
    theme_blank()
print(p)