#ifndef UTILS_VISUAL_H
#define UTILS_VISUAL_H

void set_histogram_style(TH1* h, int color, int line_width, int marker, double marker_size)
{
    h->SetLineColor(color);
    h->SetLineWidth(line_width);
    h->SetMarkerColor(color);
    h->SetMarkerStyle(marker);
    h->SetMarkerSize(marker_size);
}

#endif