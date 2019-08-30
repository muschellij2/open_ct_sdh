library(neurobase)
library(ANTsRCore)
library(extrantsr)
library(tibble)
library(dplyr)
library(tidyr)
x = list.files(
  pattern = ".jp(e|)g", path = "Patients_CT",
  recursive = TRUE, full.names = TRUE)

df = tibble(img = x,
            roi = grepl("seg", tolower(img)),
            type = basename(dirname(img)),
            id = basename(dirname(dirname(img)))) %>%
  mutate(index = basename(img),
         index = sub(".jpg", "", index),
         index = sub("_.*", "", index),
         index = as.numeric(index))
stopifnot(all(df$type %in% c("bone", "brain")))
stopifnot(!any(is.na(df$index)))
all_data = df %>% tidyr::expand(id, type, roi)

slice_range = df %>%
  group_by(id, type) %>%
  summarise(max = max(index))
slice_range = split(slice_range, list(slice_range$id, slice_range$type))
slice_range = purrr::map_df(slice_range, .f = function(x) {
  d = x %>%
    select(id, type) %>%
    distinct()
  x = tibble(id = d$id, type = d$type,
             index = 1:x$max)
  x
})

all_data = left_join(all_data, slice_range)

df = left_join(all_data, df) %>%
  arrange(id, roi, type, index) %>%
  mutate(
    outfile =
      file.path(
        "Patients_NIfTI",
        paste0(id, "_", type,
               ifelse(roi, "_ROI", "_CT"),
               ".nii.gz"))
  )

df = df %>%
  filter(!file.exists(outfile))
temp_jpg = tempfile(fileext = ".jpg")
empty_mat = matrix(0, nrow = 650, ncol = 650)
aimg = as.antsImage(empty_mat, pixeltype = "unsigned char")
antsImageWrite(aimg, temp_jpg)
aimg = antsImageRead(temp_jpg)
stopifnot(sum(aimg) == 0)

seg = df %>%
  filter(roi) %>%
  group_by(id, type) %>%
  mutate(all_missing = all(is.na(img))) %>%
  filter(!all_missing) %>%
  select(-all_missing) %>%
  ungroup()
seg = seg %>%
  mutate(img = ifelse(is.na(img), temp_jpg, img))
df = df %>%
  filter(!roi)

ss = df %>%
  ungroup() %>%
  nest(-outfile, -id, -type, -roi)
i = 1

for (i in seq(nrow(ss))) {
  print(i)
  ifile = ss[i,]
  x = ifile$data[[1]]$img
  xxx = lapply(x, antsImageRead)
  check = sapply(xxx, function(x) all(dim(x) == 650))
  stopifnot(all(check))
  xxx = lapply(xxx, as.array)

  res = abind::abind(xxx, along = 3)
  if (ifile$roi) {
    check = all(unique(res) %in% c(0, 1))
    if (!check) {
      print("not all roi 0/1")
    }
    res = res > 0
  }
  res = as.antsImage(
    res, spacing = c(0.48, 0.48, 5),
    pixeltype = ifelse(ifile$roi, "unsigned char", "float"),
    components = FALSE)
  if (!ifile$roi) {
    # This was the formula for Patients_NIfTI_minmax
    # max_value = 3071
    # min_value = -1024
    # ct = (res/255) * (max_value - min_value) + min_value
    ww = 120
    wl = 50
    min_value = wl - ww/2
    max_value = wl + ww/2
    ct = (res/255) * (max_value - min_value) + min_value
  } else {
    stopifnot(all(c(as.array(res)) %in% c(0, 1)))
    ct = res
  }

  write_nifti(ct, ifile$outfile)
}

