//
//  ATSNumericalQueue.h
//  AllTrailsShared
//
//  Created by Theodore Calmes on 5/20/13.
//  Copyright (c) 2013 AllTrails. All rights reserved.
//

#import <Foundation/Foundation.h>

/** ATSNumericalQueue is a growing queue up to a max size. It only holds float values and is mostly used for a bin for the moving average. */
@interface ATSNumericalQueue : NSObject

///-------------------------------
/// @name Queue properties
///-------------------------------

/** How many items have been pushed total. */
@property (assign, nonatomic, readonly) NSInteger pushCount;

/** Represents the maximum size the queue can reach before getting rid of its last element to accomodate for new elements. */
@property (assign, nonatomic, readonly) NSInteger size;

///-------------------------------
/// @name Initializers
///-------------------------------

/**
 @param size: The maximum number of elements the queue can hold.
 @return an ATSNumericalQueue object.
 */
- (id)initWithSize:(NSInteger)size;

///-------------------------------
/// @name Add elements
///-------------------------------

/** Add values to the queue using this method.
 @param value: a float value to append to the start of the queue.
 */
- (void)pushValue:(float)value;

///-------------------------------
/// @name Calculations using members
///-------------------------------

/** Calculates the average of the elements in the queue.
 
 Note: When the queue is smaller than the max size, it will calculate the average based on the current number of elements in the queue.
 
 @return a float value representing the average.
 */
- (float)average;

/** Calculates the standard deviation.
 @return a float representing the standard deviation.
 */
- (float)standardDeviation;

@end
