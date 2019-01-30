import sys
import datetime


def mapper_name(name):
    if name == "vg":
        return "vg"
    elif name == "vg_mitty":
        return "vg with reads simulated by Mitty"
    elif name == "seven_bridges":
        return "Seven Bridges"
    elif name == "seven_bridges_mitty":
        return "Seven Bridges with reads simulated by Mitty"
    elif name == "bwa":
        return "Linear mapper"
    elif name == "bwa_untuned":
        return "BWA-MEM untuned"
    elif name == "two_step_graph_mapper_traversemapped":
        return "Two-step approach using traversemapper"
    elif name == "two_step_graph_mapper_linearmapped":
        return "Two-step approach using linear-mapper"
    elif name == "two_step_graph_mapper_vg":
        return "Two-step approach using vg alignments"
    else:
        raise Exception("Mapper name %s not supported. Add a color for this mapper" % name)

def mapper_color(name):
    if name == "vg":
        return "#3355FF"
    elif name == "vg_mitty":
        return "#79B4FF"
    elif name == "seven_bridges":
        return "#BC35C2"
    elif name == "seven_bridges_mitty":
        return "#D797DA"
    elif name == "bwa":
        return "#B81B1B"
    elif name == "bwa_untuned":
        return "#D68B8B"
    elif name == "two_step_graph_mapper_traversemapped":
        return "#0CC1D3"
    elif name == "two_step_graph_mapper_linearmapped":
        return "#89D2D9"
    elif name == "two_step_graph_mapper_vg":
        return "#D8DE5E"
    else:
        raise Exception("Mapper name %s not supported. Add a color for this mapper" % name)

print("""
<html>
<head>
    <link rel="stylesheet" href="https://stackpath.bootstrapcdn.com/bootstrap/4.2.1/css/bootstrap.min.css" integrity="sha384-GJzZqFGwb1QTTN6wy59ffF1BuGJpLSa9DkKMp0DgiMDm4iYMj70gZWKYbI706tWS" crossorigin="anonymous">

    <style>
        .image_box{
            display: table-cell;
            padding: 10px;
        
    </style>
</head>

<body>

<div style='width: 1800px;'>
""")

images = [("All reads", ""), ("Reads with variants", "-novel"), ("Reads without variants", "-known")]

header = "Report %s" % str(datetime.datetime.now())
print("<div align='center' style='width: 1800px; font-size: 2em; font-weight: bold;'>" + header + "</div>")
directory = "./"
for title, image_postfix in images:
    image_url = "roc" + image_postfix + "-builder.png"
    width = 590
    if "novel" in image_postfix:
        width = 550

    print("""<div class='image_box'>
    <div align='center'><p style='font-size: 1.5em'><b>%s</b></p></div>
    <div style='background-size: cover; width: %dpx; height: 600px; background-image: url(%s)'></div>
    </div>
    """ % (title, width, image_url)
    )



print("<div style='font-size: 1.5em; text-align: center;'>")
for name in sys.argv[1].split(","):
    color = mapper_color(name)
    print("<span style='margin-right: 30px'><font color='" + color + "'>&#9644; " + mapper_name(name) + "</font></span>")

print("</div>")



print("""</div></body>
</html>""")


