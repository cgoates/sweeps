import numpy as np
class QuadMesh:
    def __init__(self, meshFile):
        points, quads = self.getMeshValues(meshFile)
        self.points = points
        self.quads = quads

    def __eq__(self, other):
        if not all(np.array_equal(a, b) for a, b in zip(self.points, other.points)):
            return False
        if not all(a == b for a, b in zip(self.quads, other.quads)):
            return False
        return True

    def getMeshValues(self, meshFile):
        f = open(meshFile, "r")
        lines = f.readlines()
        f.close()
        numPoints = -1
        pointsStartIndex = -1
        numQuads = -1
        quadsStartIndex = -1
        for line in range(len(lines)):
            if lines[line].startswith("POINTS"):
                pointsLine = lines[line].split()
                numPoints = int(pointsLine[1])
                pointsStartIndex = line
            elif lines[line].startswith("CELLS"):
                quadsLine = lines[line].split()
                numQuads = int(quadsLine[1])
                quadsStartIndex = line
        if (numPoints < 0) or (numQuads < 0) or (pointsStartIndex < 0) or (quadsStartIndex < 0):
            raise Exception("Given Vtk file is empty or formatted incorrectly.")
        pointLines = lines[pointsStartIndex+1:pointsStartIndex+numPoints+1]
        quadLines = lines[quadsStartIndex+1:quadsStartIndex+numQuads+1]
        points = self.getPoints(pointLines)
        quads = self.getQuads(quadLines)

        return points, quads

    def getPoints(self, lines):
        points = []
        for line in lines:
            nums = line.split()
            points.append(np.array(nums, dtype=float))
        return points

    def getQuads(self, lines):
        quads = []
        for line in lines:
            nums = line.split()
            quads.append(list(map(int, nums[1:])))
        return quads