"""
Copyright (c) 2016 Cyrill Jauner

Permission is hereby granted, free of charge, to any person obtaining a copy of this software
and associated documentation files (the "Software"), to deal in the Software without restriction,
including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense,
and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so,
subject to the following conditions:

The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO
THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL
THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT,
TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
"""
import math
from Tkinter import*
from Graphs import *


def show_graphs(G1, G2, show_numbers):
    """
    Shows the two given graphs in one frame.
    :param G1: Graph 1.
    :param G2: Graph 2.
    :param show_numbers: Whether the vertex numbers will be drawn.
    """
    master = Tk()
    master.wm_title('Graph visualization')

    # Sizes
    canvas_width = 1500
    canvas_height = 700

    w = Canvas(master,
               width=canvas_width,
               height=canvas_height)
    w.pack()

    # Draws both graphs onto the canvas.
    draw_graph(G1, show_numbers, 200, 400, 'Random graph', w)
    draw_graph(G2, show_numbers, 200, 1000,'Isomorphic copy', w)

    master.mainloop()


def draw_graph(G, show_numbers, radius, padding_left, label,  w):
    """
    Draws the given graph onto the given canvas. The vertices are circular drawn.
    :param G: The graph to draw.
    :param show_numbers: Whether the vertex numbers will be drawn.
    :param radius: The radius of the circle.
    :param padding_left: The left padding of the circle.
    :param label: The label of the graph.
    :param w: The canvas to draw onto.
    """

    nv = G.num_vertices()

    # The radius of the vertex number circle
    radius_text = radius + 40

    # Calculates the gap between two vertices on the circle.
    # There are differences between even and odd numbers.
    if nv % 2 == 0:
        step = math.pi / (nv / 2)
    else:
        step = 2*math.pi / nv

    # Writes the label of the graph.
    w.create_text(padding_left, radius-70, text=label)

    # Calculates the point on the canvas for each vertex.
    # The same calculation for the vertex labels.
    points = []
    for i in range(1,nv+1):
        y = radius*math.sin(step * i)
        x = radius*math.cos(step * i)
        x += padding_left
        y += 2 * radius
        points.append((x,y))
        w.create_oval(x-5,y-5,x+5,y+5)

        if show_numbers:
            y = radius_text * math.sin(step * i)
            x = radius_text * math.cos(step * i)
            x += 2 * radius
            y += 2 * radius
            w.create_text(x,y,text=str(i))

    # Draws each edge
    for i in range(0, len(G.edges)):
        edge = G.edges[i]
        w.create_line(points[edge[0] -1][0], points[edge[0] -1][1], points[edge[1] -1][0], points[edge[1] -1][1])


#
# Creates a random graph and an isomorphic copy.
#
g = random_graph(99,2)
g2 = g.isomorphic_copy()[0]
show_graphs(g, g2, False)