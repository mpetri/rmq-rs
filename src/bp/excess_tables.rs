pub const FWD_MIN_IDX: [u8; 256] = [
    8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 0, 8, 8, 8, 8, 8, 8, 8, 0, 8, 8, 8, 0, 2, 0, 1, 0,
    8, 8, 8, 8, 8, 8, 8, 0, 8, 8, 8, 0, 2, 0, 1, 0, 4, 4, 4, 0, 2, 0, 1, 0, 3, 3, 1, 0, 2, 0, 1, 0,
    6, 6, 6, 6, 6, 6, 6, 0, 6, 6, 6, 0, 2, 0, 1, 0, 4, 4, 4, 0, 2, 0, 1, 0, 3, 3, 1, 0, 2, 0, 1, 0,
    5, 5, 5, 5, 5, 5, 1, 0, 3, 3, 1, 0, 2, 0, 1, 0, 4, 4, 4, 0, 2, 0, 1, 0, 3, 3, 1, 0, 2, 0, 1, 0,
    7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 1, 0, 7, 7, 7, 7, 7, 7, 1, 0, 3, 3, 1, 0, 2, 0, 1, 0,
    5, 5, 5, 5, 5, 5, 1, 0, 3, 3, 1, 0, 2, 0, 1, 0, 4, 4, 4, 0, 2, 0, 1, 0, 3, 3, 1, 0, 2, 0, 1, 0,
    6, 6, 6, 6, 6, 6, 6, 0, 6, 6, 6, 0, 2, 0, 1, 0, 4, 4, 4, 0, 2, 0, 1, 0, 3, 3, 1, 0, 2, 0, 1, 0,
    5, 5, 5, 5, 5, 5, 1, 0, 3, 3, 1, 0, 2, 0, 1, 0, 4, 4, 4, 0, 2, 0, 1, 0, 3, 3, 1, 0, 2, 0, 1, 0,
];

pub const FWD_EXC: [i8; 256] = [
    -8, -6, -6, -4, -6, -4, -4, -2, -6, -4, -4, -2, -4, -2, -2, 0, -6, -4, -4, -2, -4, -2, -2, 0,
    -4, -2, -2, 0, -2, 0, 0, 2, -6, -4, -4, -2, -4, -2, -2, 0, -4, -2, -2, 0, -2, 0, 0, 2, -4, -2,
    -2, 0, -2, 0, 0, 2, -2, 0, 0, 2, 0, 2, 2, 4, -6, -4, -4, -2, -4, -2, -2, 0, -4, -2, -2, 0, -2,
    0, 0, 2, -4, -2, -2, 0, -2, 0, 0, 2, -2, 0, 0, 2, 0, 2, 2, 4, -4, -2, -2, 0, -2, 0, 0, 2, -2,
    0, 0, 2, 0, 2, 2, 4, -2, 0, 0, 2, 0, 2, 2, 4, 0, 2, 2, 4, 2, 4, 4, 6, -6, -4, -4, -2, -4, -2,
    -2, 0, -4, -2, -2, 0, -2, 0, 0, 2, -4, -2, -2, 0, -2, 0, 0, 2, -2, 0, 0, 2, 0, 2, 2, 4, -4, -2,
    -2, 0, -2, 0, 0, 2, -2, 0, 0, 2, 0, 2, 2, 4, -2, 0, 0, 2, 0, 2, 2, 4, 0, 2, 2, 4, 2, 4, 4, 6,
    -4, -2, -2, 0, -2, 0, 0, 2, -2, 0, 0, 2, 0, 2, 2, 4, -2, 0, 0, 2, 0, 2, 2, 4, 0, 2, 2, 4, 2, 4,
    4, 6, -2, 0, 0, 2, 0, 2, 2, 4, 0, 2, 2, 4, 2, 4, 4, 6, 0, 2, 2, 4, 2, 4, 4, 6, 2, 4, 4, 6, 4,
    6, 6, 8,
];

pub const FWD_MIN: [u8; 256] = [
    8, 6, 6, 4, 6, 4, 4, 2, 6, 4, 4, 2, 4, 2, 2, 0, 6, 4, 4, 2, 4, 2, 2, 0, 4, 2, 2, 0, 2, 0, 1, 0,
    6, 4, 4, 2, 4, 2, 2, 0, 4, 2, 2, 0, 2, 0, 1, 0, 4, 2, 2, 0, 2, 0, 1, 0, 3, 1, 1, 0, 2, 0, 1, 0,
    6, 4, 4, 2, 4, 2, 2, 0, 4, 2, 2, 0, 2, 0, 1, 0, 4, 2, 2, 0, 2, 0, 1, 0, 3, 1, 1, 0, 2, 0, 1, 0,
    5, 3, 3, 1, 3, 1, 1, 0, 3, 1, 1, 0, 2, 0, 1, 0, 4, 2, 2, 0, 2, 0, 1, 0, 3, 1, 1, 0, 2, 0, 1, 0,
    7, 5, 5, 3, 5, 3, 3, 1, 5, 3, 3, 1, 3, 1, 1, 0, 5, 3, 3, 1, 3, 1, 1, 0, 3, 1, 1, 0, 2, 0, 1, 0,
    5, 3, 3, 1, 3, 1, 1, 0, 3, 1, 1, 0, 2, 0, 1, 0, 4, 2, 2, 0, 2, 0, 1, 0, 3, 1, 1, 0, 2, 0, 1, 0,
    6, 4, 4, 2, 4, 2, 2, 0, 4, 2, 2, 0, 2, 0, 1, 0, 4, 2, 2, 0, 2, 0, 1, 0, 3, 1, 1, 0, 2, 0, 1, 0,
    5, 3, 3, 1, 3, 1, 1, 0, 3, 1, 1, 0, 2, 0, 1, 0, 4, 2, 2, 0, 2, 0, 1, 0, 3, 1, 1, 0, 2, 0, 1, 0,
];
