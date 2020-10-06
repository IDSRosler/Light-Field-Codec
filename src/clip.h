
#ifndef CLIP_H
#define CLIP_H

/* clip.hh */
/* Author: Pekka Astola */
/* <pekka.astola@tuni.fi>*/

template<class T>
T clip(T in, const T min, const T max) {
    if (in > max)
        return max;
    if (in < min)
        return min;
    return in;
}

#endif //CLIP_H
