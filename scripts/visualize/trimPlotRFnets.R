library(magick)

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