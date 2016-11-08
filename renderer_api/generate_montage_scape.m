function im = generate_montage_scape(rc, z, scale)
% generates a montage scape image using collection rc and section z at scale "scale"
[Wbox, bbox, url] = get_section_bounds_renderer(rc, z);
[im, v, url] = get_image_box_renderer(rc, z, Wbox, scale, 'kk');