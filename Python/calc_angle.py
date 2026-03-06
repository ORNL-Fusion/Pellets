import numpy as np

def calc_angle(entry, axis, shift):

    [x1,y1] = entry
    [x2,y2] = axis

    a = y2 - y1
    b = x2 - x1
    c = np.sqrt((x2 - x1)**2 + (y2 - y1)**2)

    print(a, b, c)

    # orig_angle = np.pi/2 - np.atan(b/a)
    # orig_angle = np.acos(a/c)
    orig_angle = np.pi/2.0 - np.atan2(b,a)
    orig_deg = orig_angle*180./np.pi
    shift_deg = orig_deg + shift
    print(shift_deg)
    shift_rad = np.pi*shift/180.0

    print(orig_angle*180./np.pi,shift_rad*180./np.pi)

    # shift_plus = [x1 - (-b*np.cos(shift_rad) + a*np.sin(shift_rad)), y2 - (-b*np.sin(shift_rad) - a*np.cos(shift_rad))]
    # shift_minus = [x1 - (-b*np.cos(shift_rad-np.pi) - a*np.sin(shift_rad-np.pi)), y2 - (b*np.sin(shift_rad-np.pi) - a*np.cos(shift_rad-np.pi))]

    # shift_plus = [x1 + (b*np.cos(shift) + a*np.sin(shift)), y2 - (-b*np.sin(shift) + a*np.cos(shift))]
    # shift_minus = [x2 - (b*np.cos(shift) + a*np.sin(shift)), y2 + (-b*np.sin(shift) + a*np.cos(shift))]

    shift_minus = [x1 + (b*np.cos(shift_rad) + a*np.sin(shift_rad)), y1 + (-b*np.sin(shift_rad) + a*np.cos(shift_rad))]
    shift_plus = [x1 + (b*np.cos(-shift_rad) + a*np.sin(-shift_rad)), y1 + (-b*np.sin(-shift_rad) + a*np.cos(-shift_rad))]

    return [orig_angle,shift_plus,shift_minus]


    