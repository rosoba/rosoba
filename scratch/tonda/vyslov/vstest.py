'''
Created on Oct 6, 2014

@author: rch
'''

import codecs

b_slova = ['byt',
           'bydlit',
           'bydlet',
           'obyvatel',
           'byt',
           'byk',
           ]


# f = codecs.open('/home/rch/Trash/text_priklad.txt', 'r', 'utf-8')


f = open('/home/rch/Trash/text_priklad.txt', 'r')

text = f.read()

words = text.split()

lower_words = []
for word in words:
    lower_words.append(word.lower())

lower_words.sort()

y_lists = dict(by=[],
               fy=[],
               ly=[],
               my=[],
               py=[],
               sy=[],
               vy=[],
               zy=[])

for word in lower_words:
    for y, y_list in y_lists.items():
        if y in word:
            y_list.append(word)

for y, y_list in y_lists.items():
    print
    print 'Vyjmenovana slova', y
    for word in y_list:
        print word,
    print
