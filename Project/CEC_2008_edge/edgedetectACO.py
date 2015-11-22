import numpy as np
import time

try:
    import matplotlib.pylab as plt
    import matplotlib.image as mpimg
except ImportError:
    pass


def acoedgedetect():
    global v
    plt.close(all)
    clear(all)
    clc
    # image loading
    filename = 'eye'
    img = rgb2gray(mpimg.imread(np.array(np.hstack((filename + '.bmp')))))
    img = np.double(img) / 255.
    [nrow, ncol] = img.shape
    # visiblity function initialization, see equation (4)
    for nMethod in range(1., 5.):
        # Four kernel functions used in the paper, see equations (7)-(10)
        # E: exponential; F: flat; G: gaussian; S:Sine; T:Turkey; W:Wave
        # fprintf('Ant Colony Optimization for Edge detection.\nPlease wait......\n');
        v = np.zeros(img.shape)
        v_norm = 0.
        for rr in range(1., nrow + 1):
            for cc in range(1., ncol + 1):
                # definition of clique
                temp1 = np.array(np.vstack((np.hstack((rr - 2., cc - 1.)), np.hstack((rr - 2., cc + 1.)),
                                            np.hstack((rr - 1., cc - 2.)), np.hstack((rr - 1., cc - 1.)),
                                            np.hstack((rr - 1., cc)), np.hstack((rr - 1., cc + 1.)),
                                            np.hstack((rr - 1., cc + 2.)), np.hstack((rr, cc - 1.)))))
                temp2 = np.array(np.vstack((np.hstack((rr + 2., cc + 1.)), np.hstack((rr + 2., cc - 1.)),
                                            np.hstack((rr + 1., cc + 2.)), np.hstack((rr + 1., cc + 1.)),
                                            np.hstack((rr + 1., cc)), np.hstack((rr + 1., cc - 1.)),
                                            np.hstack((rr + 1., cc - 2.)), np.hstack((rr, cc + 1.)))))
                temp0 = nonzero(np.logical_and(np.logical_and(np.logical_and(np.logical_and(np.logical_and(
                    np.logical_and(np.logical_and(temp1[:, 0] >= 0., temp1[:, 0] < nrow), temp1[:, 1] >= 0.),
                    temp1[:, 1] < ncol), temp2[:, 0] >= 0.), temp2[:, 0] < nrow), temp2[:, 1] >= 0.),
                    temp2[:, 1] < ncol))
                temp11 = temp1[int(temp0) - 1, :]
                temp22 = temp2[int(temp0) - 1, :]
                temp00 = np.zeros(temp11.size)
                for kk in np.arange(1, temp11.size + 1):
                    temp00[int(kk) - 1] = np.abs((
                                                 img[int(temp11[int(kk) - 1, 0]) - 1, int(temp11[int(kk) - 1, 1]) - 1] -
                                                 img[int(temp22[int(kk) - 1, 0]) - 1, int(temp22[int(kk) - 1, 1]) - 1]))

                if temp11.size == 0.:
                    v[int(rr) - 1, int(cc) - 1] = 0.
                    v_norm = v_norm + v[int(rr) - 1, int(cc) - 1]
                else:
                    lambda_val = 1.0
                    _switch_val = nMethod
                    if _switch_val == 1.:
                        # 'F'
                        temp00 *= lambda_val
                    elif _switch_val == 2.:
                        # 'Q'
                        temp00 = lambda_val * temp00 ** 2.
                    elif _switch_val == 3.:
                        # 'S'
                        temp00 = np.sin(np.pi * temp00 / 2. / lambda_val)
                    elif _switch_val == 4.:
                        # 'W'
                        temp00 = np.sin((np.pi * temp00 / lambda_val)) * np.pi * temp00 / lambda_val

            v[int(rr) - 1, int(cc) - 1] = np.sum(np.sum((temp00 ** 2.)))
            v_norm = v_norm + v[int(rr) - 1, int(cc) - 1]
        v /= v_norm
        # do normalization
        v *= 100.
        # pheromone function initialization
        p = 0.0001 * np.ones(img.shape)
        # parameter setting, see Section IV in CEC paper
        alpha = 1.
        # equation (4)
        beta = 0.1
        # equation (4)
        rho = 0.1
        # equation (11)
        phi = 0.05
        # equation (12), i.e., (9) in IEEE-CIM-06
        ant_total_num = 2. * np.round(np.sqrt(np.dot(nrow, ncol)))
        ant_pos_idx = np.zeros(ant_total_num, 2.)
        # record the location of ants
        # initialize the positions of ants
        #     rand('state', sum(clock));
        #     temp = rand(ant_total_num, 2);
        #     ant_pos_idx(:,1) = round(1 + (nrow-1) * temp(:,1)); %row index
        #     ant_pos_idx(:,2) = round(1 + (ncol-1) * temp(:,2)); %column index
        # initialize the position of ants at maximum value of visibility
        # function
        # [~,index]=max(v);
        max_value = max(v)
        index = v.index(max_value)
        ant_pos_idx[0:ncol - 1, 0] = np.arange(0, ncol)
        temp = np.arange(1., ncol + 1)
        temp1 = np.dot(np.ones(0, ncol), ncol)
        ant_pos_idx[int(ncol):, 0] = temp1 - temp
        ant_pos_idx[0:nrow, 1] = index
        ant_pos_idx[int(nrow):, 1] = index
        search_clique_mode = '8'
        # Figure 1
        # define the memory length, the positions in ant's memory are
        # non-admissible positions for the next movement
        if np.dot(nrow, ncol) == 128. * 128.:
            A = 40.
            memory_length = np.round((np.random.rand(1.) * (np.dot(1.15, A) - np.dot(0.85, A)) + np.dot(0.85, A)))
            # memory length
        elif np.dot(nrow, ncol) == 256. * 256.:
            A = 30.
            memory_length = np.round((np.random.rand(1.) * (np.dot(1.15, A) - np.dot(0.85, A)) + np.dot(0.85, A)))
            # memory length

        elif np.dot(nrow, ncol) == 512. * 512.:
            A = 20.
            memory_length = np.round((np.random.rand(1.) * (np.dot(1.15, A) - np.dot(0.85, A)) + np.dot(0.85, A)))
            # memory length

        # record the positions in ant's memory, convert 2D position-index (row, col) into
        # 1D position-index
        ant_memory = np.zeros(ant_total_num, memory_length)
        # System setup
        if np.dot(nrow, ncol) == 128. * 128.:
            total_step_num = 300.
            # the numbe of iterations?
        elif np.dot(nrow, ncol) == 256. * 256.:
            total_step_num = 900.

        elif np.dot(nrow, ncol) == 512. * 512.:
            total_step_num = 1500.

        total_iteration_num = 5.
        for iteration_idx in np.arange(1., total_iteration_num + 1):
            t = time.time()
            # record the positions where ant have reached in the last 'memory_length' iterations
            delta_p = np.zeros(nrow, ncol)
            for step_idx in np.arange(1., total_step_num + 1):
                delta_p_current = np.zeros(nrow, ncol)
                for ant_idx in np.arange(1., ant_total_num + 1):
                    ant_current_row_idx = ant_pos_idx[int(ant_idx) - 1, 0]
                    ant_current_col_idx = ant_pos_idx[int(ant_idx) - 1, 1]  # find the neighborhood of current position
                    if search_clique_mode == '4':
                        rr = ant_current_row_idx
                        cc = ant_current_col_idx
                        ant_search_range_temp = np.array(np.vstack((np.hstack((rr - 1., cc)), np.hstack((rr, cc + 1.)),
                                                                    np.hstack((rr + 1., cc)),
                                                                    np.hstack((rr, cc - 1.)))))
                    elif search_clique_mode == '8':
                        rr = ant_current_row_idx
                        cc = ant_current_col_idx
                        ant_search_range_temp = np.array(np.vstack((np.hstack((rr - 1., cc - 1.)),
                                                                    np.hstack((rr - 1., cc)),
                                                                    np.hstack((rr - 1., cc + 1.)),
                                                                    np.hstack((rr, cc - 1.)), np.hstack((rr, cc + 1.)),
                                                                    np.hstack((rr + 1., cc - 1.)),
                                                                    np.hstack((rr + 1., cc)),
                                                                    np.hstack((rr + 1., cc + 1.)))))

                    # remove the positions our of the image's range
                    temp = nonzero(np.logical_and(np.logical_and(
                        np.logical_and(ant_search_range_temp[:, 0] >= 1., ant_search_range_temp[:, 0] <= nrow),
                        ant_search_range_temp[:, 1] >= 1.), ant_search_range_temp[:, 1] <= ncol))
                    ant_search_range = ant_search_range_temp[int(temp) - 1, :]
                    # calculate the transit prob. to the neighborhood of current
                    # position
                    ant_transit_prob_v = np.zeros(ant_search_range.size, 1.)
                    ant_transit_prob_p = np.zeros(ant_search_range.size, 1.)
                    for kk in np.arange(1., (ant_search_range.size) + 1):
                        temp = np.dot(ant_search_range[int(kk) - 1, 0] - 1., ncol) + ant_search_range[int(kk) - 1, 1]
                        if length(nonzero((ant_memory[int(ant_idx) - 1, :] == temp))) == 0.:
                            # not in ant's memory
                            ant_transit_prob_v[int(kk) - 1] = v[
                                int(ant_search_range[int(kk) - 1, 0]) - 1, int(ant_search_range[int(kk) - 1, 1]) - 1]
                            ant_transit_prob_p[int(kk) - 1] = p[
                                int(ant_search_range[int(kk) - 1, 0]) - 1, int(ant_search_range[int(kk) - 1, 1]) - 1]
                        else:
                            # in ant's memory
                            ant_transit_prob_v[int(kk) - 1] = 0.
                            ant_transit_prob_p[int(kk) - 1] = 0.

                    # if all neighborhood are in memory, then the permissible search range is RE-calculated.
                    if np.logical_or(np.sum(np.sum(ant_transit_prob_v)) == 0.,
                                     np.sum(np.sum(ant_transit_prob_p)) == 0.):
                        for kk in np.arange(1., (ant_search_range.size) + 1):
                            temp = np.dot(ant_search_range[int(kk) - 1, 0] - 1., ncol) + ant_search_range[
                                int(kk) - 1, 1]
                            ant_transit_prob_v[int(kk) - 1] = v[
                                int(ant_search_range[int(kk) - 1, 0]) - 1, int(ant_search_range[int(kk) - 1, 1]) - 1]
                            ant_transit_prob_p[int(kk) - 1] = p[
                                int(ant_search_range[int(kk) - 1, 0]) - 1, int(ant_search_range[int(kk) - 1, 1]) - 1]

                    ant_transit_prob = ant_transit_prob_v ** alpha * ant_transit_prob_p ** beta / np.sum(
                        np.sum((ant_transit_prob_v ** alpha * ant_transit_prob_p ** beta)))
                    # generate a random number to determine the next position.
                    random.set_state(np.sum((100. * clock)))
                    temp = nonzero((np.cumsum(ant_transit_prob) >= np.random.rand(1.)), 1.)
                    ant_next_row_idx = ant_search_range[int(temp) - 1, 0]
                    ant_next_col_idx = ant_search_range[int(temp) - 1, 1]
                    if length(ant_next_row_idx) == 0.:
                        ant_next_row_idx = ant_current_row_idx
                        ant_next_col_idx = ant_current_col_idx

                    ant_pos_idx[int(ant_idx) - 1, 0] = ant_next_row_idx
                    ant_pos_idx[int(ant_idx) - 1, 1] = ant_next_col_idx
                    # record the delta_p_current
                    delta_p_current[
                        int(ant_pos_idx[int(ant_idx) - 1, 0]) - 1, int(ant_pos_idx[int(ant_idx) - 1, 1]) - 1] = 1.
                    # record the new position into ant's memory
                    if step_idx <= memory_length:
                        ant_memory[int(ant_idx) - 1, int(step_idx) - 1] = np.dot(ant_pos_idx[int(ant_idx) - 1, 0] - 1.,
                                                                                 ncol) + ant_pos_idx[
                                                                              int(ant_idx) - 1, 1]
                    elif step_idx > memory_length:
                        ant_memory[int(ant_idx) - 1, :] = np.roll(ant_memory[int(ant_idx) - 1, :],
                                                                  np.array(np.hstack((0., -1.))))
                        ant_memory[int(ant_idx) - 1, int(0) - 1] = np.dot(ant_pos_idx[int(ant_idx) - 1, 0] - 1., ncol) + \
                                                                   ant_pos_idx[int(ant_idx) - 1, 1]

                    # update the pheromone function (10) in IEEE-CIM-06
                    p = ((1. - rho) * p + rho * delta_p_current * v) * delta_p_current + p * np.abs(
                        (1. - delta_p_current))

                # end of ant_idx
                # update the pheromone function see equation (9) in IEEE-CIM-06
                delta_p = delta_p + (delta_p_current > 0.) > 0.
                p *= (1. - phi)
                # equation (9) in IEEE-CIM-06
            # end of step_idx

            # create the image at every iteration to see the difference
            T = func_seperate_two_class(p)
            # eq. (13)-(21), Calculate the threshold to seperate the edge map into two class
            mpimg.imsave(np.array(np.hstack((filename, nMethod, iteration_idx))),
                         np.uint8(np.abs(((p >= T) * 255. - 255.))), plt.gray(256.), 'bmp')
            timespent = time.time() - t
            print 'time spent for ' + iteration_idx + ' iteration: ' + timespent + '\n'
            # end of iteration_idx
            # generate edge map matrix
            # It uses pheromone function to determine edge?
        print 'Done!\n'
    # end of nMethod
    return


def func_seperate_two_class(I):
    # Local Variables: MAT, mu3, level, I, Threshold, mu2, N, mu, i, T, MBT, counts
    # Function calls: hist, sum, abs, cumsum, func_seperate_two_class
    #   ISODATA Compute global image threshold using iterative isodata method.
    #   LEVEL = ISODATA(I) computes a global threshold (LEVEL) that can be
    #   used to convert an intensity image to a binary image with IM2BW. LEVEL
    #   is a normalized intensity value that lies in the range [0, 1].
    #   This iterative technique for choosing a threshold was developed by Ridler and Calvard .
    #   The histogram is initially segmented into two parts using a starting threshold value such as 0 = 2B-1,
    #   half the maximum dynamic range.
    #   The sample mean (mf,0) of the gray values associated with the foreground pixels and the sample mean (mb,0)
    #   of the gray values associated with the background pixels are computed. A new threshold value 1 is now computed
    #   as the average of these two sample means. The process is repeated, based upon the new threshold,
    #   until the threshold value does not change any more.
    #
    # Reference :T.W. Ridler, S. Calvard, Picture thresholding using an iterative selection method,
    #            IEEE Trans. System, Man and Cybernetics, SMC-8 (1978) 630-632.
    # Convert all N-D arrays into a single column.  Convert to uint8 for
    # fastest histogram computation.
    I = I.flatten(1)
    # STEP 1: Compute mean intensity of image from histogram, set T=mean(I)
    [counts, N] = plt.hist(I, 256.)
    i = 1.
    mu = np.cumsum(counts)
    T[int(i) - 1] = matdiv(np.sum((N * counts)), mu[int(0) - 1])
    # STEP 2: compute Mean above T (MAT) and Mean below T (MBT) using T from
    # step 1
    mu2 = np.cumsum(counts[int((N <= T[int(i) - 1])) - 1])
    MBT = matdiv(np.sum((N[int((N <= T[int(i) - 1])) - 1] * counts[int((N <= T[int(i) - 1])) - 1])), mu2[int(0) - 1])
    mu3 = np.cumsum(counts[int((N > T[int(i) - 1])) - 1])
    MAT = matdiv(np.sum((N[int((N > T[int(i) - 1])) - 1] * counts[int((N > T[int(i) - 1])) - 1])), mu3[int(0) - 1])
    i = i + 1.
    T[int(i) - 1] = (MAT + MBT) / 2.
    # STEP 3 to n: repeat step 2 if T(i)~=T(i-1)
    Threshold = T[int(i) - 1]
    while np.abs((T[int(i) - 1] - T[int((i - 1.)) - 1])) >= 1.:
        mu2 = np.cumsum(counts[int((N <= T[int(i) - 1])) - 1])
        MBT = matdiv(np.sum((N[int((N <= T[int(i) - 1])) - 1] * counts[int((N <= T[int(i) - 1])) - 1])),
                     mu2[int(0) - 1])
        mu3 = np.cumsum(counts[int((N > T[int(i) - 1])) - 1])
        MAT = matdiv(np.sum((N[int((N > T[int(i) - 1])) - 1] * counts[int((N > T[int(i) - 1])) - 1])), mu3[int(0) - 1])
        i = i + 1.
        T[int(i) - 1] = (MAT + MBT) / 2.
        Threshold = T[int(i) - 1]

    # Normalize the threshold to the range [i, 1].
    level = Threshold
    return [level]
