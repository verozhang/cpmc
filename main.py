import numpy as np
from scipy.spatial.qhull import ConvexHull
import portion


def distance(a, b):
    return np.sqrt((a[0] - b[0]) ** 2 + (a[1] - b[1]) ** 2)
# distance


# return 0 when tangent, 1 when close, 2 when far, and -1 when overlapping
def relate_dist(a, b):
    dist = distance(a, b)
    if np.abs(dist - 1) <= 1e-6:
        return 0
    elif dist >= 2:
        return 2
    elif dist >= 1:
        return 1
    else:
        return -1
# relate_dist


# return b's relative angle to a, that is, angle between x axis and segment ab
def angle(a, b):
    return np.arctan2(b[1] - a[1], b[0] - a[0])
# angle


# calculate adj matrix for all circles
def calc_adj_matrix():
    for i in range(num):
        for j in range(num):
            adj_matrix[i][j] = relate_dist(circles[i], circles[j])
    return
# calc_adj_matrix


# fix adj matrix for one circle
def fix_adj_matrix(a):
    for i in range(num):
        adj_matrix[a][i] = relate_dist(circles[a], circles[i])
        adj_matrix[i][a] = relate_dist(circles[i], circles[a])
    return
# fix_adj_matrix


# return array for all neighbours (tangent) of a
def calc_neighbour(a):
    arr = []
    for i in range(num):
        if adj_matrix[a][i] == 0:
            arr.append(i)
    return arr
# calc_neighbour


# return array for all circles near (close to) a
def calc_near(a):
    arr = []
    for i in range(num):
        if adj_matrix[a][i] == 0 or adj_matrix[a][i] == 1:
            arr.append(i)
    return arr
# calc_near


# return if a circle is pinned (not movable w/o collision)
# a circle is pinned when it has at least 3 neighbours
# and each pair of closest neighbours has an angle smaller than \pi (180 degree)
def is_pinned(a):
    neighbour_arr = calc_neighbour(a)
    angle_arr = []
    for i in range(len(neighbour_arr)):
        neighbour = neighbour_arr[i]
        angle_arr.append(angle(circles[a], circles[neighbour]))

    # need at least 3 neighbours to pin a circle
    if len(angle_arr) < 3:
        return False

    angle_arr.sort()
    for i in range(len(angle_arr) - 1):
        if angle_arr[i+1] - angle_arr[i] >= np.pi:
            return False

    # check last pair
    if angle_arr[-1] - angle_arr[0] < np.pi:
        return False

    return True
# is_pinned


# randomly move mover around pivot.
def step_rotate(mover: int, pivot: int):
    # mover and pivot have to be tangent
    if adj_matrix[mover][pivot] != 0:
        print("Void step. Mover and pivot not tangent.")
        return
    # find available places around pivot
    nears = calc_near(pivot)
    available = portion.openclosed(-np.pi, np.pi)
    for near in nears:
        if near != mover:
            dist = distance(circles[pivot], circles[near])
            angl = angle(circles[pivot], circles[near])
            var = np.arccos(dist / 2)
            left = angl - var
            right = angl + var
            if left <= -np.pi:
                occupied = portion.open(-np.pi, right) | portion.openclosed(left + 2 * np.pi, np.pi)
            elif right > np.pi:
                occupied = portion.openclosed(left, np.pi) | portion.open(-np.pi, right - 2 * np.pi)
            else:
                occupied = portion.open(left, right)
            available -= occupied

    if available.empty:
        print("Void step. No available space next to the pivot.")
        return

    # TODO: this ratio should be changed into a parameter later
    stick_rate = 1
    # stick: move mover to a place sticking to an obstacle. Notice when only 1 circle around pivot stick is not doable.
    # TODO fix
    if np.random.rand() < stick_rate and available != portion.openclosed(-np.pi, np.pi):
        stick_angles = []
        for interval in available._intervals:
            stick_angles.append(interval.lower)
            stick_angles.append(interval.upper)
        if available.enclosure == portion.openclosed(-np.pi, np.pi):
            stick_angles.remove(-np.pi)
            stick_angles.remove(np.pi)
        stick_angle = stick_angles[np.random.randint(len(stick_angles))]
        return circles[pivot][0] + np.cos(stick_angle), circles[pivot][1] + np.sin(stick_angle)
    # TODO non-stick rotating
    else:
        pass
    return
# step_rotate


# calculate smallest bounding rectangle with given ratio from convex hull
# ratio: aspect ratio length/width of rectangle
def calc_bounding_rec(ratio: int):
    ch = ConvexHull(circles)
    cur_max_volume = 0
    for i in range(ch.nsimplex):
        # initialize current max distance on both directions
        cur_max_dist_axis = 0
        cur_dist_perps = []
        for j in range(num):
            cur_edge = ch.equations[i]
            dist_axis = - (cur_edge[0] * circles[j][0] + cur_edge[1] * circles[j][1] + cur_edge[2])
            if dist_axis > cur_max_dist_axis:
                cur_max_dist_axis = dist_axis
            dist_perp = cur_edge[1] * circles[j][0] - cur_edge[0] * circles[j][1]
            cur_dist_perps.append(dist_perp)
        # both "+1" since we are considering circle centres before, and diameter = 1
        cur_max_dist_axis += 1
        cur_max_dist_perp = max(cur_dist_perps) - min(cur_dist_perps) + 1
        # handle squares independently
        if ratio == 1:
            cur_volume = np.power(max(cur_max_dist_axis, cur_max_dist_perp), 2)
        # not a square
        elif cur_max_dist_axis > ratio * cur_max_dist_perp:
            cur_volume = np.power(cur_max_dist_axis, 2) * ratio
        elif cur_max_dist_perp > ratio * cur_max_dist_axis:
            cur_volume = np.power(cur_max_dist_perp, 2) * ratio
        else:
            cur_volume = np.power(min(cur_max_dist_axis, cur_max_dist_perp), 2) * ratio
        if cur_volume > cur_max_volume:
            cur_max_volume = cur_volume
    return cur_max_volume
# calc_bounding_rec


# one Metropolis Monte Carlo loop
def mc_met(k: int):
    mover = np.random.randint(num)
    if is_pinned(mover):
        mover = np.random.randint(num)
    mover_old_parameter = (circles[mover][0], circles[mover][1])
    old_energy = calc_bounding_rec(1)
    print("Old energy: ", old_energy)
    pivots = []
    for i in range(num):
        if adj_matrix[i][mover] == 0:
            pivots.append(i)
    if not pivots:
        return
    j = np.random.randint(len(pivots))
    pivot = pivots[j]

    mover_new_parameter = step_rotate(mover, pivot)
    circles[mover] = mover_new_parameter
    new_energy = calc_bounding_rec(1)
    print("New energy: ", new_energy)

    if new_energy <= old_energy:
        print("Taking step with energy falling.")
        calc_adj_matrix()
        return
    else:
        threshold = np.power(np.e, (new_energy - old_energy) * -k)
        accept = np.random.rand()
        if accept < threshold:
            print("Taking step with energy raising.")
            calc_adj_matrix()
        else:
            print("Refusing step with energy raising.")
            circles[mover] = mover_old_parameter
        return
# mc_met


num = 0
circles = np.ndarray((1,))
adj_matrix = np.ndarray((1,))


if __name__ == "__main__":
    # num = int(input("Input num of circles."))
    num = 9
    adj_matrix = np.zeros((num, num))
    circles = np.zeros((num, 2))

    circles[0][0] = 0
    circles[0][1] = 0
    circles[1][0] = 0
    circles[1][1] = 1
    circles[2][0] = 0
    circles[2][1] = 2
    circles[3][0] = 0
    circles[3][1] = 3
    circles[4][0] = 0
    circles[4][1] = 4
    circles[5][0] = 0
    circles[5][1] = 5
    circles[6][0] = 0
    circles[6][1] = 6
    circles[7][0] = 0
    circles[7][1] = 7
    circles[8][0] = 1
    circles[8][1] = 0
    '''
    circles[0][0] = 0
    circles[0][1] = 1
    circles[1][0] = 1
    circles[1][1] = 0
    circles[2][0] = 1
    circles[2][1] = 1
    circles[3][0] = 1
    circles[3][1] = 2
    circles[4][0] = 2
    circles[4][1] = 1
    '''
    cur_ch = ConvexHull(circles)
    calc_adj_matrix()
    for i in range(100):
        mc_met(1)

    pass
