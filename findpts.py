
# Fusion360API Python script

import traceback
import adsk.fusion
import adsk.core
import random

POINT_COUNT = 20

def run(context):
    ui = adsk.core.UserInterface.cast(None)
    try:
        app: adsk.core.Application = adsk.core.Application.get()
        ui = app.userInterface
        
        points = '[ '

        msg: str = 'Select Face'
        selFilter: str = 'Faces'
        face1: adsk.core.Selection = selectEnt(msg, selFilter)

        points = createPointsOnFace(
            points,
            face1.entity,
            POINT_COUNT
        )

        points += ']'
        ui.messageBox(str(len(points)))
        ui.messageBox(points[:2000])

        
    

    except:
        if ui:
            ui.messageBox('Failed:\n{}'.format(traceback.format_exc()))


def createPointsOnFace(
    points: str,
    face: adsk.fusion.BRepFace,
    count: int) -> adsk.fusion.Sketch:

    # Evaluator
    eva: adsk.core.SurfaceEvaluator = face.evaluator

    # Parameter Range
    prmRange: adsk.core.BoundingBox2D = eva.parametricRange()
    prmMin: adsk.core.Point2D = prmRange.minPoint
    prmMax: adsk.core.Point2D = prmRange.maxPoint

    # Parameters On Face
    prmsOnFace = []
    while len(prmsOnFace) < count + 1:

        prmLst = [
            (random.uniform(prmMin.x, prmMax.x), random.uniform(prmMin.y, prmMax.y))
            for _ in range(count)
        ]

        prmLst = list(set(prmLst) - set(prmsOnFace))
        for x, y in prmLst:
            prm: adsk.core.Point2D = adsk.core.Point2D.create(x, y)
            if eva.isParameterOnFace(prm):
                prmsOnFace.append((x, y))

    # Points on Face
    prms = [adsk.core.Point2D.create(x, y) for x, y in prmsOnFace]
    _, pnts = eva.getPointsAtParameters(prms)


    for pnt in pnts:
        points +=(str(int(pnt.x)/100) + ', ' + str(int(pnt.y)/100) + ', ' + str(int(pnt.z)/100) + '; ')
    return points




def selectEnt(
        msg: str,
        filterStr: str) -> adsk.core.Selection:

    try:
        app = adsk.core.Application.get()
        ui = app.userInterface
        sel = ui.selectEntity(msg, filterStr)
        return sel
    except:
        return None