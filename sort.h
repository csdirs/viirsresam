// column break points
short SORT_BREAK_POINTS[11] = {
	5, 87, 170, 358, 567, 720, 850, 997, 1120, 1275, 1600,
};

// relative row that the pixel comes from
short SORT_FIRST[16][11] = {
	{0,   0,  0,  0,  0,  0,  0,  0,  0,  0, 0},
	{0,   0,  0,  0,  0,  0,  0,  0,  0,  0, 0},
	{0,   0,  0,  0,  0,  0,  0,  0,  0,  0, 0},
	{0,   0,  0,  0,  0,  0,  0,  0,  0,  0, 0},
	{0,   0,  0,  0,  0,  0,  0,  0,  0,  0, 0},
	{0,   0,  0,  0,  0,  0,  0,  0,  0,  0, 0},
	{0,   0,  0,  0,  0,  0,  0,  0,  0,  0, 0},
	{0,   0,  0,  0,  0,  0,  0,  0,  0,  0, 0},
	{+8, +8, +8,  0,  0,  0,  0,  0,  0,  0, 0},
	{+8, -1, -1, +7,  0,  0,  0,  0,  0,  0, 0},
	{-2, +7, +7, -1, +6,  0,  0,  0,  0,  0, 0},
	{+7, -2, -2, +6, -1, +5,  0,  0,  0,  0, 0},
	{-3, +6, +6, -2, +5, -1, +4,  0,  0,  0, 0},
	{+6, -3, -3, +5, -2, +4, -1, +3,  0,  0, 0},
	{-4, +5, +5, -3, +4, -2, +3, -1, +2,  0, 0},
	{+5, -4, -4, +4, -3, +3, -2, +2, -1, +1, 0},
};

// relative row that the pixel comes from
short SORT_MID[16][11] = {
	{-5,  +4, +4, -4, +3, -3, +2, -2, +1, -1, 0},
	{+4,  -5, -5, +3, -4, +2, -3, +1, -2,  0, 0},
	{-6,  +3, +3, -5, +2, -4, +1, -3,  0,  0, 0},
	{+3,  -6, -6, +2, -5, +1, -4,  0,  0,  0, 0},
	{-7,  +2, +2, -6, +1, -5,  0,  0,  0,  0, 0},
	{+11, -7, -7, +1, -6,  0,  0,  0,  0,  0, 0},
	{+1,  +1, +1, -7,  0,  0,  0,  0,  0,  0, 0},
	{-9,  +9, -8,  0,  0,  0,  0,  0,  0,  0, 0},
	{+9,  -9, +8,  0,  0,  0,  0,  0,  0,  0, 0},
	{-1,  -1, -1, +7,  0,  0,  0,  0,  0,  0, 0},
	{-11, +7, +7, -1, +6,  0,  0,  0,  0,  0, 0},
	{+7,  -2, -2, +6, -1, +5,  0,  0,  0,  0, 0},
	{-3,  +6, +6, -2, +5, -1, +4,  0,  0,  0, 0},
	{+6,  -3, -3, +5, -2, +4, -1, +3,  0,  0, 0},
	{-4,  +5, +5, -3, +4, -2, +3, -1, +2,  0, 0},
	{+5,  -4, -4, +4, -3, +3, -2, +2, -1, +1, 0},
};

// relative row that the pixel comes from
short SORT_LAST[16][11] = {
	{-5, +4, +4, -4, +3, -3, +2, -2, +1, -1, 0},
	{+4, -5, -5, +3, -4, +2, -3, +1, -2,  0, 0},
	{-6, +3, +3, -5, +2, -4, +1, -3,  0,  0, 0},
	{+3, -6, -6, +2, -5, +1, -4,  0,  0,  0, 0},
	{-7, +2, +2, -6, +1, -5,  0,  0,  0,  0, 0},
	{+2, -7, -7, +1, -6,  0,  0,  0,  0,  0, 0},
	{-8, +1, +1, -7,  0,  0,  0,  0,  0,  0, 0},
	{-8, -8, -8,  0,  0,  0,  0,  0,  0,  0, 0},
	{0,   0,  0,  0,  0,  0,  0,  0,  0,  0, 0},
	{0,   0,  0,  0,  0,  0,  0,  0,  0,  0, 0},
	{0,   0,  0,  0,  0,  0,  0,  0,  0,  0, 0},
	{0,   0,  0,  0,  0,  0,  0,  0,  0,  0, 0},
	{0,   0,  0,  0,  0,  0,  0,  0,  0,  0, 0},
	{0,   0,  0,  0,  0,  0,  0,  0,  0,  0, 0},
	{0,   0,  0,  0,  0,  0,  0,  0,  0,  0, 0},
	{0,   0,  0,  0,  0,  0,  0,  0,  0,  0, 0},
};
