'''
Created on Jan 15, 2014

@author: rch
'''

class Point(object):

    def __init__(self):
        self.x = 0.0
        self.y = 0.0

    def move(self, dx, dy):
        self.x += dx
        self.y += dy

    def __str__(self):
        return 'x: %g, y: %g\n' % (self.x, self.y)

if __name__ == '__main__':
    p = Point()
    print 'point', p
    p.move(0.6, 4)
    print 'moved point, p', p
