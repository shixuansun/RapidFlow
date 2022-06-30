//
// Created by sunsx on 08/06/21.
//

#ifndef RAPIDMATCH_STREAMING_CONFIG_H
#define RAPIDMATCH_STREAMING_CONFIG_H

#define ENABLE_PERFORMANCE_COUNTERS

/**
 * Both modes does not list the results.
 * 0: Count the results.
 * 1: Enumerate the results.
 */
#define EXECUTION_MODE 0

/**
 * Measure the time cost for a batch of updates
 */
// #define MEASURE_UPDATE_COST
// #define MEASURE_BATCH_SIZE (128)

extern bool g_exit;
#endif //RAPIDMATCH_STREAMING_CONFIG_H
