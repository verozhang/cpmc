import numpy as np
from scipy.spatial.qhull import ConvexHull
import matplotlib.pyplot as plt
import portion


def distance(a, b):
    return np.sqrt((a[0] - b[0]) ** 2 + (a[1] - b[1]) ** 2)
# distance


# return b's relative angle to a, that is, angle between x axis and segment ab
def angle(a, b):
    return np.arctan2(b[1] - a[1], b[0] - a[0])
# angle


# return 0 when tangent, 1 when close, 2 when far, and -1 when overlapping
def calc_tangency(a, b):
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


# calculate adj matrix for all circles
def calc_adj_matrix():
    for i in range(num):
        for j in range(num):
            adj_matrix[i][j] = calc_tangency(circles[i], circles[j])
    return
# calc_adj_matrix


# fix adj matrix for one circle
def fix_adj_matrix(a):
    for i in range(num):
        adj_matrix[a][i] = calc_tangency(circles[a], circles[i])
        adj_matrix[i][a] = calc_tangency(circles[i], circles[a])
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

    # need at least 3 neighbours to pin a circle
    if len(neighbour_arr) < 3:
        return False

    angle_arr = []
    for i in range(len(neighbour_arr)):
        neighbour = neighbour_arr[i]
        angle_arr.append(angle(circles[a], circles[neighbour]))

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
def step(mover: int, pivot: int, stick_rate: float):
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
        return circles[mover]

    # stick: move mover to a place sticking to an obstacle. Notice when only 1 circle around pivot stick is not doable.
    if np.random.rand() < stick_rate:
        stick_angles = []
        for interval in available._intervals:
            stick_angles.append(interval.lower)
            stick_angles.append(interval.upper)
        if available.enclosure == portion.openclosed(-np.pi, np.pi):
            stick_angles.remove(-np.pi)
            stick_angles.remove(np.pi)
        if len(stick_angles) == 0:
            print("Not doable.")
            return circles[mover]
        stick_angle = stick_angles[np.random.randint(len(stick_angles))]
        return circles[pivot][0] + np.cos(stick_angle), circles[pivot][1] + np.sin(stick_angle)

    # non-stick: uniformly pick angle from all available angle intervals
    else:
        # available is 1 interval
        if available.atomic:
            stick_angle = available.lower + np.random.rand() * (available.upper - available.lower)
            return circles[pivot][0] + np.cos(stick_angle), circles[pivot][1] + np.sin(stick_angle)
        # available is union of 2 intervals
        else:
            stick_angle = available.lower + np.random.rand() * \
                          (available._intervals[0].upper - available._intervals[0].lower +
                           available._intervals[1].upper - available._intervals[1].lower)
            if stick_angle > available._intervals[0].upper:
                stick_angle += (available._intervals[1].lower - available._intervals[0].upper)
            return circles[pivot][0] + np.cos(stick_angle), circles[pivot][1] + np.sin(stick_angle)
            pass
# step_rotate


# calculate smallest bounding rectangle with given ratio from convex hull
# shape: aspect ratio length/width of rectangle
def calc_bounding_rec(shape: float):
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
        if shape == 1:
            cur_volume = np.power(max(cur_max_dist_axis, cur_max_dist_perp), 2)
        # not a square
        elif cur_max_dist_axis > shape * cur_max_dist_perp:
            cur_volume = np.power(cur_max_dist_axis, 2) * shape
        elif cur_max_dist_perp > shape * cur_max_dist_axis:
            cur_volume = np.power(cur_max_dist_perp, 2) * shape
        else:
            cur_volume = np.power(min(cur_max_dist_axis, cur_max_dist_perp), 2) * shape
        if cur_volume > cur_max_volume:
            cur_max_volume = cur_volume
    return cur_max_volume
# calc_bounding_rec


# one Metropolis Monte Carlo loop
# acc_parameter: parameter used to determine whether to accept an energy-raising step. Check Boltzmann distribution.
def mc_met(acc: float, shape: float, stick_rate: float):
    mover = np.random.randint(num)
    if is_pinned(mover):
        mover = np.random.randint(num)
    mover_old_parameter = (circles[mover][0], circles[mover][1])
    print("Old location: ", mover_old_parameter)
    old_energy = calc_bounding_rec(shape)
    print("Old energy: ", old_energy)

    pivots = []
    for i in range(num):
        if adj_matrix[i][mover] == 0:
            pivots.append(i)

    if not pivots or len(pivots) == 1:  # TODO split this to detect lone pair move
        ch = ConvexHull(circles)
        for j in range(len(ch.vertices)):
            pivots.append(j)
    # TODO lone pair move
    if len(pivots) == 1:
        ch = ConvexHull(circles)

    k = np.random.randint(len(pivots))
    pivot = pivots[k]

    mover_new_parameter = step(mover, pivot, stick_rate)
    circles[mover] = mover_new_parameter
    print("New location: ", mover_new_parameter)
    new_energy = calc_bounding_rec(shape)
    print("New energy: ", new_energy)

    if new_energy - old_energy <= 1e-6:
        print("Taking step with energy falling.")
        calc_adj_matrix()
        return
    else:
        threshold = np.power(np.e, (new_energy - old_energy) * -acc)
        accept = np.random.rand()
        if accept < threshold:
            print("Taking step with energy raising.")
            calc_adj_matrix()
        else:
            print("Refusing step with energy raising.")
            circles[mover] = mover_old_parameter
        return
# mc_met

# TODO shrink all components to make them touch

num = 0
circles = np.ndarray((1,))
adj_matrix = np.ndarray((1,))


if __name__ == "__main__":
    # num = int(input("Input num of circles."))
    num = 50
    adj_matrix = np.zeros((num, num))
    circles = np.zeros((num, 2))

    for i in range(num - 1):
        circles[i] = (0, i)
    circles[num - 1] = (1, 0)
    '''
    for i in range(2):
        for j in range(100):
            circles[i*100 + j] = (i, j)
    '''
    calc_adj_matrix()
    for i in range(1000000):
        print("Step ", i)
        mc_met(5, 3, 0.8)

        if i % 100 == 0:
            # plotting
            fig, ax = plt.subplots()
            fig.set_size_inches(10, 40)
            plt.xlim(-10, 10)
            plt.ylim(-10, 60)
            plt.grid(linestyle='--')
            ax.set_aspect(1)

            for i in range(num):
                cur_circle = plt.Circle(circles[i], 0.5)
                ax.add_artist(cur_circle)
            plt.show()
