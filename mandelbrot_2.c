#include <stdlib.h>
#include <stdio.h>

struct Image {
    unsigned nb_lig;
    unsigned nb_col;
    unsigned *pixels;
};

void calculer (struct Image * im, unsigned nb_iter, double x_min, double x_max, double y_min, double y_max) {
    double pasx = (x_max - x_min) / im -> nb_col;
    double pasy = (y_max - y_min) / im -> nb_lig;
    unsigned l,c;
#pragma omp parallel for private (c)
    for (l = 0; l < im->nb_lig; l++) {
        for (c = 0; c < im->nb_col; c++) {
            double zx = 0.0, zy = 0.0, new_zx;
            double cx = x_min + c*pasx, cy = y_min + l*pasy;
            unsigned n = 0;
            for(n=0;  (zx*zx + zy*zy < 4.0 ) && ( n < nb_iter ); n++ ) {
                new_zx = zx*zx - zy*zy + cx;
                zy = 2.0*zx*zy + cy;
                zx = new_zx;
            }
            if(n == nb_iter) n = 0;
            im->pixels[l*im->nb_col + c] = n;
        }
    }
}

void draw_image(struct Image *im) {
    const char charset[] = ".,c8M@jawrpogOQEPGJ";
    unsigned l,c;
    for (l = 0; l < im->nb_lig; l++) {
        for (c = 0; c < im->nb_col; c++) {
            unsigned n = im->pixels[l*im->nb_col + c];
            char p = n > 0 ? charset[n % (sizeof(charset)-1)] : ' ';
            putchar(p);
            if(c+1 == im->nb_col) puts("");
        }
    }
    puts("");
}

int main(void) {
    struct Image im;
    im.nb_lig = 40;
    im.nb_col = 80;
    im.pixels = malloc(sizeof *im.pixels * im.nb_lig*im.nb_col);
    unsigned nb_iter = 256;
    calculer(&im, nb_iter, -2.5, 1.5, -2.0, 2.0);
    draw_image(&im);
    return 0;
}
