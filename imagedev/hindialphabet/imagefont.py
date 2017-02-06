# coding=utf-8
import PIL
from PIL import ImageFont
from PIL import Image
from PIL import ImageDraw
import sys

# font = ImageFont.truetype("Arial-Bold.ttf",14)
#font = ImageFont.truetype("Arial.ttf",14)

for i in xrange(0x0900,0x097F):
    #print i
    #print 'U+'+str(i)
    #'''
    img=Image.new("RGBA", (250,250),(0,0,0))
    draw = ImageDraw.Draw(img)

    unicode_font = ImageFont.truetype("Shree714.ttc", 200)
    #draw.text ( (10,10), unicode_text, font=unicode_font, fill=font_color )


    draw.text((5, 5),unichr(i),(255,255,255),font=unicode_font)
    draw = ImageDraw.Draw(img)
    img.save("letters/"+unichr(i)+".png")
    img.save("letters/"+"{:02x}".format(i)+".png")

    #sys.exit()
    #'''
