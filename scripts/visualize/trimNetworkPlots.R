library(magick)

## RF-Net figures

# (1) trim out white-space bordering figure, (2) add back a thinner border,
# (3) save image
l2r_7ret_fig <- image_read_pdf("figures/l2r_7ret.pdf")
image_write(
  image_border(image_trim(l2r_7ret_fig),"white","100x100"),
  "figures/l2r_7ret_trim.pdf",format="pdf",density=600)

r2l_6ret_fig <- image_read_pdf("figures/r2l_6ret.pdf")
image_write(
  image_border(image_trim(r2l_6ret_fig),"white","100x100"),
  "figures/r2l_6ret_trim.pdf",format="pdf",density=600)

l2r_7ret_97g_fig <- image_read_pdf("figures/l2r_7ret_97g.pdf")
image_write(
  image_border(image_trim(l2r_7ret_97g_fig),"white","100x100"),
  "figures/l2r_7ret_97g_trim.pdf",format="pdf",density=600)

r2l_7ret_97g_fig <- image_read_pdf("figures/r2l_7ret_97g.pdf")
image_write(
  image_border(image_trim(r2l_7ret_97g_fig),"white","100x100"),
  "figures/r2l_7ret_97g_trim.pdf",format="pdf",density=600)

###############################################################################

## NetRAX figures

l2r_15g_no3cyc_fig <- image_read_pdf("figures/netrax/multistart-set1c-L2R-15genes-no3cyc.pdf")
image_write(
  image_border(image_trim(l2r_15g_no3cyc_fig),"white","100x100"),
  "figures/netrax/multistart-set1c-L2R-15genes-no3cyc_trim.pdf",format="pdf",density=600)

r2l_15g_no3cyc_fig <- image_read_pdf("figures/netrax/multistart-set1c-R2L-15genes-no3cyc.pdf")
image_write(
  image_border(image_trim(r2l_15g_no3cyc_fig),"white","100x100"),
  "figures/netrax/multistart-set1c-R2L-15genes-no3cyc_trim.pdf",format="pdf",density=600)

l2r_1b_15g_3cyc_fig <- image_read_pdf("figures/netrax/multistart-set1b-L2R-15genes-with3cyc.pdf")
image_write(
  image_border(image_trim(l2r_1b_15g_3cyc_fig),"white","100x100"),
  "figures/netrax/multistart-set1b-L2R-15genes-with3cyc_trim.pdf",format="pdf",density=600)

l2r_1b_15g_no3cyc_fig <- image_read_pdf("figures/netrax/multistart-set1b-L2R-15genes-no3cyc.pdf")
image_write(
  image_border(image_trim(l2r_1b_15g_no3cyc_fig),"white","100x100"),
  "figures/netrax/multistart-set1b-L2R-15genes-no3cyc_trim.pdf",format="pdf",density=600)
