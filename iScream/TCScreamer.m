//
//  TCScreamer.m
//  iScream
//
//  Created by Theodore Calmes on 5/2/14.
//  Copyright (c) 2014 theo. All rights reserved.
//

#import "TCScreamer.h"
#import "ATSNumericalQueue.h"

@import AVFoundation;

@interface TCScreamer ()

@property (strong) CMMotionManager *manager;
@property BOOL playing;
@property (strong) AVAudioPlayer *player;
@property (strong) ATSNumericalQueue *queue;

@end

@implementation TCScreamer

- (instancetype)init
{
    self = [super init];
    if (self) {
        self.manager = [[CMMotionManager alloc] init];
        self.playing = NO;
        self.player = [[AVAudioPlayer alloc] initWithContentsOfURL:[[NSBundle mainBundle] URLForResource:@"falling" withExtension:@"wav"] error:nil];
        self.queue = [[ATSNumericalQueue alloc] initWithSize:3];
        [self.manager startDeviceMotionUpdatesToQueue:[NSOperationQueue new] withHandler:^(CMDeviceMotion *motion, NSError *error) {
            CMAcceleration a = motion.userAcceleration;
            double mag = sqrt(a.x * a.x + a.y * a.y + a.z * a.z);
            [self.queue pushValue:mag];
            if ([self.queue average] >= 0.3 && !self.playing) {
                [[NSOperationQueue mainQueue] addOperationWithBlock:^{
                    [self.player play];
                }];
                self.playing = YES;
            } else if (self.playing && [self.queue average] < 0.3) {
                [[NSOperationQueue mainQueue] addOperationWithBlock:^{
                    [self.player stop];
                }];
                self.playing = NO;
            }
        }];
    }
    return self;
}

@end
