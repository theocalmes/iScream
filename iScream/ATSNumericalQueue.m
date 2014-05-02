//
//  ATSNumericalQueue.m
//  AllTrailsShared
//
//  Created by Theodore Calmes on 5/20/13.
//  Copyright (c) 2013 AllTrails. All rights reserved.
//

#import "ATSNumericalQueue.h"
#include <Accelerate/Accelerate.h>
#import "ATSMath.h"

@interface ATSNumericalQueue ()
@property (assign, nonatomic) float *values;
@end

@implementation ATSNumericalQueue

- (id)initWithSize:(NSInteger)size
{
    NSParameterAssert(size != 0);
    
    self = [super init];
    if (self) {
        _size = size;
        _values = (float *)calloc(_size, sizeof(float));
        _pushCount = 0;
    }
    return self;
}

- (void)pushValue:(float)value
{
    float *newValues = (float *)calloc(_size, sizeof(float));
    for (int i = 1; i < _size; i++) {
        newValues[i] = _values[i - 1];
    }
    newValues[0] = value;

    for (NSInteger i = 0; i < _size; i++) {
        _values[i] = newValues[i];
    }

    free(newValues);

    _pushCount++;
}

- (float)average
{
    return ATSAverage(_values, MIN(_size, _pushCount));
}

- (float)standardDeviation
{
    return ATSStandardDeviation(_values, MIN(_size, _pushCount));
}

- (void)dealloc
{
    free(_values);
}

@end
