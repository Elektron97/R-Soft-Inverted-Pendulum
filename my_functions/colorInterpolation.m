function color_faded = colorInterpolation(init_color, end_color, n, i)
    alpha = i/n;
    color_faded = alpha*end_color + (1 - alpha)*init_color;
end