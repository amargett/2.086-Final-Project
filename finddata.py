#Author-
#Description-

import adsk.core, adsk.fusion, adsk.cam, traceback

def run(context):
    ui: adsk.core.UserInterface = None
    try:
        _app : adsk.core.Application = adsk.core.Application.get()
        _ui = _app.userInterface
        des: adsk.fusion.Design = _app.activeProduct
        root: adsk.fusion.Component = des.rootComponent

        selFiltter: str = 'Edges, SketchLines'
        front = _ui.selectEntity('Select front line', selFiltter)
        bowline = _ui.selectEntity('Select bow line', selFiltter)
        sternline = _ui.selectEntity('Select stern line', selFiltter)
        back = _ui.selectEntity('Select back line', selFiltter)
        side1 = _ui.selectEntity('Select side 1, near chine line', selFiltter)
        side2 = _ui.selectEntity('Select side 2', selFiltter)
        bottomline = _ui.selectEntity('Select bottom line', selFiltter)
        topline = _ui.selectEntity('Select top line', selFiltter)
        chineline = _ui.selectEntity('Select chine line', selFiltter)


                # Minimum Distance

        measMgr = _app.measureManager
        LOA =  adsk.core.MeasureResults = measMgr.measureMinimumDistance(front.entity, back.entity)
        LB = adsk.core.MeasureResults = measMgr.measureMinimumDistance(front.entity, bowline.entity)
        LS = adsk.core.MeasureResults = measMgr.measureMinimumDistance(sternline.entity, back.entity)
        BD = adsk.core.MeasureResults = measMgr.measureMinimumDistance(side1.entity, side2.entity)
        DD = adsk.core.MeasureResults = measMgr.measureMinimumDistance(topline.entity, bottomline.entity)
        BC = adsk.core.MeasureResults = measMgr.measureMinimumDistance(chineline.entity, side2.entity)
        DC = adsk.core.MeasureResults = measMgr.measureMinimumDistance(chineline.entity, bottomline.entity)



        # result
        _ui.messageBox('LOA: ' + str(LOA.value) + 'LB: ' + str(LB.value) + 'LS: ' + str(LS.value) + 'BD: ' + str(BD.value) + 'DD: ' + str(DD.value)
            + 'BC: ' + str(BC.value) + 'DC: ' + str(DC.value))



    except:
        if _ui:
            _ui.messageBox('Failed:\n{}'.format(traceback.format_exc()))


