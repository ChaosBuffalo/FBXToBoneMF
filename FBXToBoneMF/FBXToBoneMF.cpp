#include <iostream>
#include <fbxsdk.h>
#undef snprintf;
#include "json.hpp"
#include <vector>
#include <map>
#include <set>
#include <fstream>
#include <iomanip>
#include <sstream>
#include <ios>
#define _USE_MATH_DEFINES
#include <math.h>


using nlohmann::json;


int numTabs = 0;


void PrintTabs() {
	for (int i = 0; i < numTabs; i++)
		printf("\t");
}

std::string GetRotationOrder(EFbxRotationOrder rotationOrder) {
	switch (rotationOrder) {
	case EFbxRotationOrder::eOrderXYZ:
		return "xyz";
	case EFbxRotationOrder::eOrderXZY:
		return "xzy";
	case EFbxRotationOrder::eOrderYXZ:
		return "yxz";
	case EFbxRotationOrder::eOrderYZX:
		return "yzx";
	case EFbxRotationOrder::eOrderZXY:
		return "zxy";
	case EFbxRotationOrder::eOrderSphericXYZ:
		return "spheric_zxy";
	default:
		return "unknown";
	}
}

/**
 * Return a string-based representation based on the attribute type.
 */
FbxString GetAttributeTypeName(FbxNodeAttribute::EType type) {
	switch (type) {
	case FbxNodeAttribute::eUnknown: return "unidentified";
	case FbxNodeAttribute::eNull: return "null";
	case FbxNodeAttribute::eMarker: return "marker";
	case FbxNodeAttribute::eSkeleton: return "skeleton";
	case FbxNodeAttribute::eMesh: return "mesh";
	case FbxNodeAttribute::eNurbs: return "nurbs";
	case FbxNodeAttribute::ePatch: return "patch";
	case FbxNodeAttribute::eCamera: return "camera";
	case FbxNodeAttribute::eCameraStereo: return "stereo";
	case FbxNodeAttribute::eCameraSwitcher: return "camera switcher";
	case FbxNodeAttribute::eLight: return "light";
	case FbxNodeAttribute::eOpticalReference: return "optical reference";
	case FbxNodeAttribute::eOpticalMarker: return "marker";
	case FbxNodeAttribute::eNurbsCurve: return "nurbs curve";
	case FbxNodeAttribute::eTrimNurbsSurface: return "trim nurbs surface";
	case FbxNodeAttribute::eBoundary: return "boundary";
	case FbxNodeAttribute::eNurbsSurface: return "nurbs surface";
	case FbxNodeAttribute::eShape: return "shape";
	case FbxNodeAttribute::eLODGroup: return "lodgroup";
	case FbxNodeAttribute::eSubDiv: return "subdiv";
	default: return "unknown";
	}
}

struct BoneMFVertex {
	double x;
	double y;
	double z;
	double u;
	double v;
	double nX;
	double nY;
	double nZ;
	std::vector<std::pair<std::string, double>>* boneWeights;
	bool normalsAdded;

	bool isEqual(BoneMFVertex other) {
		return nX == other.nX && nY == other.nY && nZ == other.nZ && u == other.u
			&& v == other.v && x == other.x && y == other.y && z == other.z;
	}

	BoneMFVertex(const BoneMFVertex& rhs) {
		x = rhs.x;
		y = rhs.y;
		z = rhs.z;
		nX = rhs.nX;
		nY = rhs.nY;
		nZ = rhs.nZ;
		u = rhs.u;
		v = rhs.v;
		normalsAdded = rhs.normalsAdded;
		boneWeights = new std::vector<std::pair<std::string, double>>(*rhs.boneWeights);
	}

	BoneMFVertex& operator=(const BoneMFVertex& rhs) {
		if (this != &rhs) {
			x = rhs.x;
			y = rhs.y;
			z = rhs.z;
			nX = rhs.nX;
			nY = rhs.nY;
			nZ = rhs.nZ;
			u = rhs.u;
			v = rhs.v;
			normalsAdded = rhs.normalsAdded;
			delete boneWeights;
			boneWeights = new std::vector<std::pair<std::string, double>>(*rhs.boneWeights);
		}
		return *this;
	}

	BoneMFVertex(double x, double y, double z)
		: x(x), y(y), z(z)
	{
		boneWeights = new std::vector<std::pair<std::string, double>>();
		normalsAdded = false;
	}

	void addBoneWeight(std::string boneName, double weight) {
		boneWeights->push_back(std::make_pair(boneName, weight));
	}

	~BoneMFVertex()
	{
		delete boneWeights;
	}
};




struct BoneMFTriangle {
	int vertexIndices[3];
};

struct BoneMFNodeFrame {
	FbxVector4 trans;
	FbxVector4 rot;
	FbxVector4 scale;
};

struct BoneMFAnimationChannel {
	std::vector<BoneMFNodeFrame>* frames;

	BoneMFAnimationChannel()
	{
		frames = new std::vector<BoneMFNodeFrame>();

	}

	void addFrame(BoneMFNodeFrame frame) {
		frames->push_back(frame);
	}

	BoneMFAnimationChannel(const BoneMFAnimationChannel& rhs) {
		frames = new std::vector<BoneMFNodeFrame>(*rhs.frames);
	}

	BoneMFAnimationChannel& operator=(const BoneMFAnimationChannel& rhs) {

		if (this != &rhs) {
			delete frames;
			frames = new std::vector<BoneMFNodeFrame>(*rhs.frames);
		}
		return *this;
	}


	~BoneMFAnimationChannel()
	{
		delete frames;
	}
};



struct Mesh {

	std::vector<BoneMFTriangle> triangles;
};


BoneMFVertex createVertex(double x, double y, double z) {

	return BoneMFVertex(x, y, z);;
}

BoneMFTriangle createTriangle(int ind1, int ind2, int ind3) {
	BoneMFTriangle triangle;
	triangle.vertexIndices[0] = ind1;
	triangle.vertexIndices[1] = ind2;
	triangle.vertexIndices[2] = ind3;
	return triangle;
}

/**
 * Print an attribute.
 */

json vertexToJson(const BoneMFVertex& vertex) {
	json vertexJson;
	vertexJson["x"] = vertex.x;
	vertexJson["y"] = vertex.y;
	vertexJson["z"] = vertex.z;
	vertexJson["u"] = vertex.u;
	vertexJson["v"] = vertex.v;
	vertexJson["nX"] = vertex.nX;
	vertexJson["nY"] = vertex.nY;
	vertexJson["nZ"] = vertex.nZ;
	json boneWeights;
	for (auto& pair : *vertex.boneWeights) {
		json weight;
		weight["boneName"] = pair.first;
		weight["weight"] = pair.second;
		boneWeights.push_back(weight);
	}
	vertexJson["boneWeights"] = boneWeights;
	return vertexJson;
}

json matrixToJson(const FbxAMatrix& mat) {
	json matrix;
	for (int i = 0; i < 4; i++) {
		json rowJson;
		FbxVector4 row = mat.GetRow(i);
		rowJson.push_back(row[0]);
		rowJson.push_back(row[1]);
		rowJson.push_back(row[2]);
		rowJson.push_back(row[3]);
		matrix.push_back(rowJson);
	}
	return matrix;
}

json triangleToJson(const BoneMFTriangle& triangle) {
	json triangleJson;
	triangleJson["indices"] = triangle.vertexIndices;
	return triangleJson;
}

json meshToJson(
	std::vector<BoneMFVertex>& vertices,
	std::vector<BoneMFTriangle>& triangles) {
	json jsonVertices;
	json jsonTriangles;
	for (BoneMFVertex& vert : vertices) {
		jsonVertices.push_back(vertexToJson(vert));
	}
	for (BoneMFTriangle& tri : triangles) {
		jsonTriangles.push_back(tri.vertexIndices[0]);
		jsonTriangles.push_back(tri.vertexIndices[1]);
		jsonTriangles.push_back(tri.vertexIndices[2]);
	}
	json mesh;
	mesh["vertices"] = jsonVertices;
	mesh["triangles"] = jsonTriangles;

	return mesh;


}

void ProcessMesh(FbxMesh* mesh, std::vector<BoneMFVertex>& vertices,
	std::vector<BoneMFTriangle>& triangles) {
	int vertexCount = mesh->GetControlPointsCount();
	int polygonCount = mesh->GetPolygonCount();

	for (int i = 0; i < vertexCount; i++) {
		FbxVector4 point = mesh->GetControlPointAt(i);

		vertices.push_back(createVertex(point.mData[0], point.mData[1], point.mData[2]));
	}

	for (int i = 0; i < polygonCount; i++) {
		int polySize = mesh->GetPolygonSize(i);
		BoneMFTriangle triangle = createTriangle(mesh->GetPolygonVertex(i, 0),
			mesh->GetPolygonVertex(i, 1), mesh->GetPolygonVertex(i, 2));
		triangles.push_back(triangle);
	}

	// process bone weights
	int deformerCount = mesh->GetDeformerCount();
	for (int i = 0; i < deformerCount; i++) {
		FbxDeformer* deformer = mesh->GetDeformer(i);
		if (deformer->GetDeformerType() == FbxDeformer::eSkin) {
			FbxSkin* skin = (FbxSkin*)deformer;
			for (int iCluster = 0; iCluster < skin->GetClusterCount(); iCluster++) {
				FbxCluster* cluster = skin->GetCluster(iCluster);
				FbxNode* boneNode = cluster->GetLink();
				if (boneNode == NULL) {
					printf("Cluster has no bone, skipping \n");
					continue;
				}
				else {
					printf("Cluster %d has bone: %s \n", iCluster, boneNode->GetName());
				}
				switch (cluster->GetLinkMode()) {
				case FbxCluster::eNormalize: {
					printf("Normalize Link \n");
					break;
				}
				case FbxCluster::eAdditive: {
					printf("Additive Link \n");
					break;
				}
				case FbxCluster::eTotalOne: {
					printf("Total Link \n");
					break;
				}
				default: {
					printf("Invalid link mode \n");
				}
				}

				int iNumVertsAttached = cluster->GetControlPointIndicesCount();
				printf("%d vertices attached \n", iNumVertsAttached);

				double* pdWeights = cluster->GetControlPointWeights();
				int* piVertIndex = cluster->GetControlPointIndices();

				for (int iV = 0; iV < iNumVertsAttached; ++iV) {
					BoneMFVertex& vertex = vertices.at(piVertIndex[iV]);
					vertex.addBoneWeight(boneNode->GetName(), pdWeights[iV]);
					printf("%i{%f}, ", piVertIndex[iV], pdWeights[iV]);
				}
				printf("\n\n");
			}
		}
	}

	std::map<int, std::vector<int>> vertexDuplicates;

	for (int i = 0; i < polygonCount; i++) {
		int polySize = mesh->GetPolygonSize(i);
		for (int c = 0; c < polySize; c++) {
			FbxVector4 normal;
			mesh->GetPolygonVertexNormal(i, c, normal);
			int polyVertIndex = mesh->GetPolygonVertex(i, c);
			printf("<vertex index='%d' x='%f' y='%f' z='%f', w='%f' \n",
				polyVertIndex, normal.mData[0], normal.mData[1], normal.mData[2], normal.mData[3]);
			auto& vertex = vertices.at(polyVertIndex);
			if (vertex.normalsAdded) {
				printf("Normals arent equal %d, %d, %d \n", i, c, polyVertIndex);
				BoneMFVertex newVertex = BoneMFVertex(vertex);
				newVertex.nX = normal.mData[0];
				newVertex.nY = normal.mData[1];
				newVertex.nZ = normal.mData[2];
				vertices.push_back(newVertex);
				int newIndex = vertices.size() - 1;
				auto& triangle = triangles.at(i);
				triangle.vertexIndices[c] = newIndex;
			}
			else {
				vertex.nX = normal.mData[0];
				vertex.nY = normal.mData[1];
				vertex.nZ = normal.mData[2];
				vertex.normalsAdded = true;
			}

		}
	}



	// Get UV sets
	FbxStringList lUVSetNameList;
	mesh->GetUVSetNames(lUVSetNameList);

	// Get the first uv set, we won't support different coordinates in BoneTown
	const char* lUVSetName = lUVSetNameList.GetStringAt(0);
	printf("UV set: '%s' \n", lUVSetName);
	const FbxGeometryElementUV* lUVElement = mesh->GetElementUV(lUVSetName);
	const bool lUseIndex = lUVElement->GetReferenceMode() != FbxGeometryElement::eDirect;
	const int lIndexCount = (lUseIndex) ? lUVElement->GetIndexArray().GetCount() : 0;

	if (lUVElement->GetMappingMode() == FbxGeometryElement::eByControlPoint)
	{
		for (int lPolyIndex = 0; lPolyIndex < polygonCount; ++lPolyIndex)
		{
			// build the max index array that we need to pass into MakePoly
			const int lPolySize = mesh->GetPolygonSize(lPolyIndex);
			for (int lVertIndex = 0; lVertIndex < lPolySize; ++lVertIndex)
			{
				FbxVector2 lUVValue;

				//get the index of the current vertex in control points array
				int lPolyVertIndex = mesh->GetPolygonVertex(lPolyIndex, lVertIndex);
				int vertIndex = triangles.at(lPolyIndex).vertexIndices[lVertIndex];
				//the UV index depends on the reference mode
				int lUVIndex = lUseIndex ? lUVElement->GetIndexArray().GetAt(lPolyVertIndex) : lPolyVertIndex;

				lUVValue = lUVElement->GetDirectArray().GetAt(lUVIndex);

				auto& vertex = vertices.at(vertIndex);
				vertex.u = lUVValue.mData[0];
				vertex.v = lUVValue.mData[1];
			}
		}
	}
	else if (lUVElement->GetMappingMode() == FbxGeometryElement::eByPolygonVertex)
	{
		int lPolyIndexCounter = 0;
		for (int lPolyIndex = 0; lPolyIndex < polygonCount; ++lPolyIndex)
		{
			// build the max index array that we need to pass into MakePoly
			const int lPolySize = mesh->GetPolygonSize(lPolyIndex);
			for (int lVertIndex = 0; lVertIndex < lPolySize; ++lVertIndex)
			{
				if (lPolyIndexCounter < lIndexCount)
				{
					int lPolyVertIndex = mesh->GetPolygonVertex(lPolyIndex, lVertIndex);
					FbxVector2 lUVValue;
					int vertIndex = triangles.at(lPolyIndex).vertexIndices[lVertIndex];
					//the UV index depends on the reference mode
					int lUVIndex = lUseIndex ? lUVElement->GetIndexArray().GetAt(lPolyIndexCounter) : lPolyIndexCounter;

					lUVValue = lUVElement->GetDirectArray().GetAt(lUVIndex);

					auto& vertex = vertices.at(vertIndex);
					vertex.u = lUVValue.mData[0];
					vertex.v = lUVValue.mData[1];
					lPolyIndexCounter++;
				}
			}
		}
	}

	//std::set<int> duplicateSet;


	//int i = 0;
	//std::vector<std::pair<int, int>> duplicates;
	//for (auto& vertex : vertices) {
	//	int c = 0;
	//	for (auto& other : vertices) {
	//		if (vertex.isEqual(other) && c != i && duplicateSet.find(i) == duplicateSet.end()) {
	//			printf("Found 2 equal verts %d, %d \n ", c, i);
	//			duplicates.push_back(std::make_pair(c, i));
	//			duplicateSet.insert(c);
	//		}
	//		c++;
	//	}
	//	i++;
	//}

	//std::vector<BoneMFVertex> finalVertices;

	//int triCount = 0;
	//for (auto& triangle : triangles) {
	//	printf("triangle: %d, ind1: %d, ind2: %d, ind3: %d \n", triCount, triangle.vertexIndices[0],
	//		triangle.vertexIndices[1], triangle.vertexIndices[2]);
	//	triCount++;
	//}

	//// Fix triangles to only reference the original

	//for (auto& duplicate : duplicates) {
	//	triCount = 0;
	//	for (BoneMFTriangle& triangle : triangles) {
	//		for (int triIndex = 0; triIndex < 3; triIndex++) {
	//			if (triangle.vertexIndices[triIndex] == duplicate.first) {
	//				printf("Swapping triangle %d : %d from %d to %d \n", triCount, triIndex, triangle.vertexIndices[triIndex], duplicate.second);
	//				triangle.vertexIndices[triIndex] = duplicate.second;
	//			}
	//		}
	//		triCount++;
	//	}

	//}

	//triCount = 0;
	//for (auto& triangle : triangles) {
	//	printf("triangle: %d, ind1: %d, ind2: %d, ind3: %d \n", triCount, triangle.vertexIndices[0],
	//		triangle.vertexIndices[1], triangle.vertexIndices[2]);
	//	triCount++;
	//}

	//std::map<int, int> oldToNew;
	//// Remove duplicates from vertices
	//int count = 0;
	//for (auto& vertex : vertices) {
	//	if (duplicateSet.find(count) != duplicateSet.end()) {
	//		printf("Found dupliate for %d \n", count);
	//	}
	//	else {
	//		finalVertices.push_back(vertex);
	//		oldToNew.insert(std::make_pair(count, finalVertices.size() - 1));
	//		printf("Vertex mapping %d -> %d  (%d) \n", count, oldToNew.at(count), finalVertices.size() - 1);
	//	}
	//	count++;
	//}

	//// Fix triangles to reference new indices
	//triCount = 0;
	//for (BoneMFTriangle& triangle : triangles) {
	//	for (int triIndex = 0; triIndex < 3; triIndex++) {
	//		int originalVertex = triangle.vertexIndices[triIndex];
	//		if (oldToNew.find(originalVertex) != oldToNew.end()) {
	//			printf("Swapping triangle %d:%d from %d to %d \n", triCount, triIndex, originalVertex, oldToNew.at(originalVertex));
	//			triangle.vertexIndices[triIndex] = oldToNew.at(originalVertex);
	//		}
	//	}
	//	triCount++;
	//}

	int i = 0;
	for (auto& vertex : vertices) {
		printf("<vertex index='%d' x='%f' y='%f' z='%f', u='%f', v='%f', nX='%f', nY='%f', nZ='%f' \n", i,
			vertex.x, vertex.y, vertex.z, vertex.u, vertex.v, vertex.nX, vertex.nY, vertex.nZ);
		i++;
	}

	int triCount = 0;
	for (auto& triangle : triangles) {
		printf("triangle: %d, ind1: %d, ind2: %d, ind3: %d \n", triCount, triangle.vertexIndices[0],
			triangle.vertexIndices[1], triangle.vertexIndices[2]);
		triCount++;
	}


}

json attributeToJson(FbxNodeAttribute* pAttribute) {
	json attributeJson;
	if (!pAttribute) {
		return attributeJson;
	}
	FbxString typeName = GetAttributeTypeName(pAttribute->GetAttributeType());
	attributeJson["type"] = typeName;
	switch (pAttribute->GetAttributeType()) {
	case FbxNodeAttribute::eMesh: {

		FbxMesh* pMesh = (FbxMesh*)pAttribute;
		std::vector<BoneMFVertex> vertices;
		std::vector<BoneMFTriangle> triangles;
		ProcessMesh(pMesh, vertices, triangles);
		attributeJson["mesh"] = meshToJson(vertices, triangles);
		break;
	}
	}
	FbxString attrName = pAttribute->GetName();
	PrintTabs();
	// Note: to retrieve the character array of a FbxString, use its Buffer() method.
	printf("<attribute type='%s' name='%s'/>\n", typeName.Buffer(), attrName.Buffer());
	return attributeJson;
}


double degreesToRadians(double degrees) {
	return (degrees * M_PI) / 180.0;
}

void printMatrix(FbxAMatrix matrix, std::string name) {
	printf("'%s' matrix \n", name.c_str());
	for (int i = 0; i < 4; i++) {
		FbxVector4 vec = matrix.GetColumn(i);
		printf("Row: %d: %f, %f, %f, %f \n", i, vec[0], vec[1], vec[2], vec[3]);
	}
}

void printVector(FbxVector4 vec, std::string name) {
	printf("'%s': %f, %f, %f, %f \n", name.c_str(), vec[0], vec[1], vec[2], vec[3]);
}

std::map<std::string, FbxAMatrix> visited;

FbxAMatrix CalculateGlobalTransform(FbxNode* pNode)
{
	FbxAMatrix lTranslationM, lScalingM, lScalingPivotM, lScalingOffsetM, lRotationOffsetM, lRotationPivotM, \
		lPreRotationM, lRotationM, lPostRotationM, lTransform;
	FbxAMatrix lParentGX, lGlobalT, lGlobalRS;
	if (!pNode)
	{
		lTransform.SetIdentity();
		return lTransform;
	}
	if (visited.find(pNode->GetName()) != visited.end()) {
		return visited.at(pNode->GetName());
	}

	// Construct translation matrix
	FbxVector4 lTranslation = pNode->LclTranslation.Get();
	lTranslationM.SetT(lTranslation);
	// Construct rotation matrices
	FbxVector4 lRotation = pNode->LclRotation.Get();
	FbxVector4 lPreRotation = pNode->PreRotation.Get();
	FbxVector4 lPostRotation = pNode->PostRotation.Get();
	lRotationM.SetR(lRotation);
	lPreRotationM.SetR(lPreRotation);
	lPostRotationM.SetR(lPostRotation);
	// Construct scaling matrix
	FbxVector4 lScaling = pNode->LclScaling.Get();
	lScalingM.SetS(lScaling);
	// Construct offset and pivot matrices
	FbxVector4 lScalingOffset = pNode->ScalingOffset.Get();
	FbxVector4 lScalingPivot = pNode->ScalingPivot.Get();
	FbxVector4 lRotationOffset = pNode->RotationOffset.Get();
	FbxVector4 lRotationPivot = pNode->RotationPivot.Get();
	lScalingOffsetM.SetT(lScalingOffset);
	lScalingPivotM.SetT(lScalingPivot);
	lRotationOffsetM.SetT(lRotationOffset);
	lRotationPivotM.SetT(lRotationPivot);

	printf("Starting Node Transform Calculation: '%s' \n", pNode->GetName());
	printf("Rotation order: %s \n", GetRotationOrder(pNode->RotationOrder).c_str());
	printMatrix(lTranslationM, "translationM");
	printVector(lRotation, "vRotation");
	printMatrix(lRotationM, "rotationM");
	printVector(lPostRotation, "vPostRotation");
	printMatrix(lPostRotationM, "postRotationM");
	printVector(lPreRotation, "vPreRotation");
	printMatrix(lPreRotationM, "preRotationM");
	printMatrix(lScalingM, "scalingM");
	printMatrix(lScalingOffsetM, "scalingOffsetM");
	printMatrix(lScalingPivotM, "scalingPivotM");
	printMatrix(lRotationOffsetM, "rotationOffsetM");
	printMatrix(lRotationPivotM, "rotationPivotM");
	// Calculate the global transform matrix of the parent node
	FbxNode* lParentNode = pNode->GetParent();
	if (lParentNode)
	{
		lParentGX = CalculateGlobalTransform(lParentNode);
	}
	else
	{
		lParentGX.SetIdentity();
	}

	printMatrix(lParentGX, "parent Global");
	//Construct Global Rotation
	FbxAMatrix lLRM, lParentGRM;
	FbxVector4 lParentGR = lParentGX.GetR();
	lParentGRM.SetR(lParentGR);
	printMatrix(lParentGRM, "parent global rotation");
	lLRM = lPreRotationM * lRotationM * lPostRotationM;
	printMatrix(lLRM, "local rot mat");
	//Construct Global Shear*Scaling
	//FBX SDK does not support shear, to patch this, we use:
	//Shear*Scaling = RotationMatrix.Inverse * TranslationMatrix.Inverse * WholeTranformMatrix
	FbxAMatrix lLSM, lParentGSM, lParentGRSM, lParentTM;
	FbxVector4 lParentGT = lParentGX.GetT();
	lParentTM.SetT(lParentGT);
	printMatrix(lParentTM, "parent global trans");

	lParentGRSM = lParentTM.Inverse() * lParentGX;
	lParentGSM = lParentGRM.Inverse() * lParentGRSM;

	printMatrix(lParentGRSM, "parent global rot scale");
	printMatrix(lParentGSM, "parent global scale");
	lLSM = lScalingM;

	//Do not consider translation now
	FbxTransform::EInheritType lInheritType = pNode->InheritType.Get();
	if (lInheritType == FbxTransform::eInheritRrSs)
	{
		lGlobalRS = lParentGRM * lLRM * lParentGSM * lLSM;
	}
	else if (lInheritType == FbxTransform::eInheritRSrs)
	{
		lGlobalRS = lParentGRM * lParentGSM * lLRM * lLSM;
	}
	else if (lInheritType == FbxTransform::eInheritRrs)
	{
		FbxAMatrix lParentLSM;
		FbxVector4 lParentLS = lParentNode->LclScaling.Get();
		lParentLSM.SetS(lParentLS);
		printMatrix(lParentLSM, "parent local scaling");
		FbxAMatrix lParentGSM_noLocal = lParentGSM * lParentLSM.Inverse();
		printMatrix(lParentGSM_noLocal, "parent global scalign no local");
		lGlobalRS = lParentGRM * lLRM * lParentGSM_noLocal * lLSM;
	}
	else
	{
		FBXSDK_printf("error, unknown inherit type! \n");
	}

	printMatrix(lGlobalRS, "global rot scale");

	printMatrix(lRotationPivotM.Inverse(), "rot pivot inverse");
	printMatrix(lScalingPivotM.Inverse(), "scale pivot inverse");

	FbxAMatrix tempRotate = lRotationOffsetM * lRotationPivotM * lPreRotationM * lRotationM * lPostRotationM * lRotationPivotM.Inverse();
	FbxAMatrix tempScaling = lScalingOffsetM * lScalingPivotM * lScalingM * lScalingPivotM.Inverse();
	printMatrix(tempRotate, "rotation multiplier");
	printMatrix(tempScaling, "scaling multiplier");
	// Construct translation matrix
	// Calculate the local transform matrix
	lTransform = lTranslationM * lRotationOffsetM * lRotationPivotM * lPreRotationM * lRotationM * lPostRotationM * lRotationPivotM.Inverse()\
		* lScalingOffsetM * lScalingPivotM * lScalingM * lScalingPivotM.Inverse();
	printMatrix(lTransform, "local transform");
	FbxVector4 lLocalTWithAllPivotAndOffsetInfo = lTransform.GetT();
	printVector(lLocalTWithAllPivotAndOffsetInfo, "pivot/offset translate");
	// Calculate global translation vector according to: 
	// GlobalTranslation = ParentGlobalTransform * LocalTranslationWithPivotAndOffsetInfo
	FbxVector4 lGlobalTranslation = lParentGX.MultT(lLocalTWithAllPivotAndOffsetInfo);
	printVector(lGlobalTranslation, "global translation");
	lGlobalT.SetT(lGlobalTranslation);
	printMatrix(lGlobalT, "global translation m");
	printMatrix(lGlobalRS, "global rot scale");
	//Construct the whole global transform
	lTransform = lGlobalT * lGlobalRS;
	visited.insert(std::make_pair(pNode->GetName(), lTransform));
	return lTransform;
}

json parseNodeAnimations(FbxNode* node) {
	json animationJson;
	return animationJson;
}
/**
 * Print a node, its attributes, and all its children recursively.
 */
json nodeToJson(FbxNode* pNode, FbxScene* scene) {
	json nodeJson;
	PrintTabs();
	const char* nodeName = pNode->GetName();
	nodeJson["name"] = nodeName;
	FbxVector4 translation = pNode->LclTranslation.Get();
	nodeJson["translation"] = { translation[0], translation[1], translation[2], translation[3] };
	FbxVector4  rotation = pNode->LclRotation.Get();
	nodeJson["rotation"] = { degreesToRadians(rotation[0]), degreesToRadians(rotation[1]),
		degreesToRadians(rotation[2]), rotation[3] };
	FbxVector4  preRotation = pNode->PreRotation.Get();
	nodeJson["preRotation"] = { degreesToRadians(preRotation[0]), degreesToRadians(preRotation[1]),
		degreesToRadians(preRotation[2]), preRotation[3] };
	FbxVector4  postRotation = pNode->PostRotation.Get();
	nodeJson["postRotation"] = { degreesToRadians(postRotation[0]), degreesToRadians(postRotation[1]),
		degreesToRadians(postRotation[2]), postRotation[3] };
	FbxVector4  scaling = pNode->LclScaling.Get();
	nodeJson["scaling"] = { scaling[0], scaling[1], scaling[2], scaling[3] };
	FbxVector4  scalingOffset = pNode->ScalingOffset.Get();
	nodeJson["scalingOffset"] = { scalingOffset[0], scalingOffset[1], scalingOffset[2], scalingOffset[3] };
	FbxVector4  scalingPivot = pNode->ScalingPivot.Get();
	nodeJson["scalingPivot"] = { scalingPivot[0], scalingPivot[1], scalingPivot[2], scalingPivot[3] };
	FbxVector4  rotationOffset = pNode->RotationOffset.Get();
	nodeJson["rotationOffset"] = { rotationOffset[0], rotationOffset[1], rotationOffset[2], rotationOffset[3] };
	FbxVector4  rotationPivot = pNode->RotationPivot.Get();
	nodeJson["rotationPivot"] = { rotationPivot[0], rotationPivot[1], rotationPivot[2], rotationPivot[3] };


	auto inheritType = pNode->InheritType.Get();
	switch (inheritType) {
	case FbxTransform::EInheritType::eInheritRrSs:
	{
		nodeJson["inheritType"] = "RrSs";
		break;
	}
	case FbxTransform::EInheritType::eInheritRSrs:
	{
		nodeJson["inheritType"] = "RSrs";
		break;
	}
	case FbxTransform::EInheritType::eInheritRrs:
	{
		nodeJson["inheritType"] = "Rrs";
		break;
	}
	default:
	{
		nodeJson["inheritType"] = "unknown";
	}
	}
	FbxAMatrix lGlobal, lLocal;
	lGlobal = pNode->EvaluateGlobalTransform();
	//nodeJson["global"] = matrixToJson(lGlobal);
	lLocal = pNode->EvaluateLocalTransform();
	//nodeJson["local"] = matrixToJson(lLocal);

	//FbxAMatrix calculatedGlobal = CalculateGlobalTransform(pNode);

	//printf("Global Transform: %s \n", nodeName);
	//for (int i = 0; i < 4; i++) {
	//	FbxVector4 vec = calculatedGlobal.GetColumn(i);
	//	printf("Row: %d: %f, %f, %f, %f \n", i, vec[0], vec[1], vec[2], vec[3]);
	//}

	//printf("Evaluate Global Transform: %s \n", nodeName);
	//for (int i = 0; i < 4; i++) {
	//	FbxVector4 vec = lGlobal.GetColumn(i);
	//	printf("Row: %d: %f, %f, %f, %f \n", i, vec[0], vec[1], vec[2], vec[3]);
	//}

	// Print the contents of the node.
	printf("<node name='%s' translation='(%f, %f, %f)' rotation='(%f, %f, %f)' scaling='(%f, %f, %f)'>\n",
		nodeName,
		translation[0], translation[1], translation[2],
		rotation[0], rotation[1], rotation[2],
		scaling[0], scaling[1], scaling[2]
	);
	numTabs++;

	// Print the node's attributes.
	json attributesJson;
	for (int i = 0; i < pNode->GetNodeAttributeCount(); i++)
	{
		attributesJson.push_back(attributeToJson(pNode->GetNodeAttributeByIndex(i)));
	}

	nodeJson["attributes"] = attributesJson;


	json childrenNodes;
	// Recursively jsonify the children
	for (int j = 0; j < pNode->GetChildCount(); j++) {
		childrenNodes.push_back(nodeToJson(pNode->GetChild(j), scene));
	}

	nodeJson["children"] = childrenNodes;

	numTabs--;
	PrintTabs();

	printf("</node>\n");

	return nodeJson;
}

void parseAnimLayer(FbxAnimLayer* animLayer, FbxNode* node,
	std::map<std::string, BoneMFAnimationChannel>& channels) {

	BoneMFAnimationChannel channel;

	channels.insert(std::make_pair(node->GetName(), channel));

}

FbxTime getStartTime(FbxScene* pScene, FbxAnimStack* animStack) {
	FbxTakeInfo* lCurrentTakeInfo = pScene->GetTakeInfo(animStack->GetName());

	FbxTime mStart;
	if (lCurrentTakeInfo)
	{
		mStart = lCurrentTakeInfo->mLocalTimeSpan.GetStart();
	}
	else
	{
		// Take the time line value
		FbxTimeSpan lTimeLineTimeSpan;
		pScene->GetGlobalSettings().GetTimelineDefaultTimeSpan(lTimeLineTimeSpan);

		mStart = lTimeLineTimeSpan.GetStart();
	}

	return mStart;
}

FbxTime getAnimationDuration(FbxScene* pScene, FbxAnimStack* animStack) {


	FbxTakeInfo* lCurrentTakeInfo = pScene->GetTakeInfo(animStack->GetName());

	FbxTime mStart, mStop;
	if (lCurrentTakeInfo)
	{
		mStart = lCurrentTakeInfo->mLocalTimeSpan.GetStart();
		mStop = lCurrentTakeInfo->mLocalTimeSpan.GetStop();
	}
	else
	{
		// Take the time line value
		FbxTimeSpan lTimeLineTimeSpan;
		pScene->GetGlobalSettings().GetTimelineDefaultTimeSpan(lTimeLineTimeSpan);

		mStart = lTimeLineTimeSpan.GetStart();
		mStop = lTimeLineTimeSpan.GetStop();
	}
	printf("Start time: %s stopTime: %s \n", mStart.GetTimeString(), mStop.GetTimeString());
	FbxTime animationDuration = mStop - mStart;
	return animationDuration;
}

void parseNodeForFrames(FbxNode* node, FbxAnimEvaluator* evaluator, FbxLongLong frameCount,
	std::map<std::string, BoneMFAnimationChannel>& channels, FbxTime frameStart) {
	printf("Parsing node for animations: %s \n", node->GetName());
	BoneMFAnimationChannel channel;
	for (FbxLongLong i = 0; i < frameCount; i++) {
		FbxTime time;
		time.SetFrame(i);
		time += frameStart;
		BoneMFNodeFrame frame;
		FbxVector4 lclRot = evaluator->GetNodeLocalRotation(node, time);
		FbxVector4 lclScale = evaluator->GetNodeLocalScaling(node, time);
		FbxVector4 lclTrans = evaluator->GetNodeLocalTranslation(node, time);
		frame.rot = lclRot;
		frame.scale = lclScale;
		frame.trans = lclTrans;
		channel.addFrame(frame);
		//printVector(lclRot, "rot");
		//printVector(lclScale, "scale");
		//printVector(lclTrans, "trans");
	}
	channels.insert(std::make_pair(node->GetName(), channel));
	for (int j = 0; j < node->GetChildCount(); j++) {
		parseNodeForFrames(node->GetChild(j), evaluator, frameCount, channels, frameStart);
	}
}

json nodeFrameToJson(BoneMFNodeFrame& frame) {
	json frameJson;
	frameJson["rotation"] = { degreesToRadians(frame.rot[0]), degreesToRadians(frame.rot[1]),
		degreesToRadians(frame.rot[2]), frame.rot[3] };
	frameJson["translation"] = { frame.trans[0], frame.trans[1],
		frame.trans[2], frame.trans[3] };
	frameJson["scale"] = { frame.scale[0], frame.scale[1], frame.scale[2], frame.scale[3] };
	return frameJson;

}

json channelToJson(BoneMFAnimationChannel& channel) {
	json channelJson;
	for (auto& frame : *channel.frames) {
		channelJson.push_back(nodeFrameToJson(frame));
	}
	return channelJson;
}


json parseAnimStack(FbxScene* pScene, FbxAnimStack* animStack, FbxNode* rootNode) {
	FbxAnimEvaluator* evaluator = pScene->GetAnimationEvaluator();
	pScene->SetCurrentAnimationStack(animStack);
	int nbAnimLayers = animStack->GetMemberCount<FbxAnimLayer>();
	printf("Found %d animation layers \n", nbAnimLayers);
	FbxTime time = getAnimationDuration(pScene, animStack);
	FbxTime frameStart = getStartTime(pScene, animStack);
	FbxLongLong frameCount = time.GetFrameCount();
	double frameRate = time.GetFrameRate(time.GetGlobalTimeMode());
	printf("Time for anim: %s \n", time.GetTimeString());
	printf("Frame rate: %f \n", time.GetFrameRate(time.GetGlobalTimeMode()));
	printf("Frame count: %d \n", time.GetFrameCount());
	std::map<std::string, BoneMFAnimationChannel> channels;
	parseNodeForFrames(rootNode, evaluator, frameCount, channels, frameStart);
	json animJson;
	animJson["frameRate"] = frameRate;
	animJson["frameCount"] = frameCount;
	json channelsJson;
	for (std::map<std::string, BoneMFAnimationChannel>::iterator it = channels.begin(); it != channels.end(); ++it) {
		channelsJson[it->first] = channelToJson(it->second);
	}
	animJson["channels"] = channelsJson;
	return animJson;

}

/**
 * Main function - loads the hard-coded fbx file,
 * and prints its contents in an xml format to stdout.
 */
int main(int argc, char** argv) {

	// Change the following filename to a suitable filename value.

	char* lFilename = argv[1];
	char* newName = argv[2];
	char* fileFormat = argv[3];
	std::string fileF(fileFormat);
	bool doJson = true;
	if (fileF == "cbor") {
		doJson = false;
	}
	else if (fileF == "json") {
		doJson = true;
	}
	else {
		printf("Only json and cbor file formats are supported. Exiting.");
		return 0;
	}
	printf("Filename %s \n", lFilename);
	// Initialize the SDK manager. This object handles all our memory management.
	FbxManager* lSdkManager = FbxManager::Create();

	// Create the IO settings object.
	FbxIOSettings* ios = FbxIOSettings::Create(lSdkManager, IOSROOT);
	lSdkManager->SetIOSettings(ios);

	// Create an importer using the SDK manager.
	FbxImporter* lImporter = FbxImporter::Create(lSdkManager, "");

	// Use the first argument as the filename for the importer.
	if (!lImporter->Initialize(lFilename, -1, lSdkManager->GetIOSettings())) {
		printf("Call to FbxImporter::Initialize() failed.\n");
		printf("Error returned: %s\n\n", lImporter->GetStatus().GetErrorString());
		exit(-1);
	}

	// Create a new scene so that it can be populated by the imported file.
	FbxScene* lScene = FbxScene::Create(lSdkManager, "myScene");

	// Import the contents of the file into the scene.
	lImporter->Import(lScene);


	// The file is imported; so get rid of the importer.
	lImporter->Destroy();


	// Triangulate meshes
	FbxGeometryConverter converter(lSdkManager);
	converter.Triangulate(lScene, true);
	converter.SplitMeshesPerMaterial(lScene, true);


	FbxNode* lRootNode = lScene->GetRootNode();
	json modelFile;
	json nodes;
	if (lRootNode) {
		for (int i = 0; i < lRootNode->GetChildCount(); i++) {
			FbxNode* childNode = lRootNode->GetChild(i);
			json node = nodeToJson(childNode, lScene);
			nodes.push_back(node);
		}
	}
	modelFile["nodes"] = nodes;

	int animations = lScene->GetSrcObjectCount(FbxCriteria::ObjectType(FbxAnimStack::ClassId));
	printf("Number animations %d \n", animations);

	json animationJson;
	for (int i = 0; i < animations; i++)
	{
		FbxAnimStack* lAnimStack = lScene->GetSrcObject<FbxAnimStack>(i);
		FbxString lOutputString = "Animation Stack Name: ";
		lOutputString += lAnimStack->GetName();
		lOutputString += "\n\n";
		printf(lOutputString);
		animationJson[lAnimStack->GetName()] = parseAnimStack(lScene, lAnimStack, lRootNode);
	}
	modelFile["animations"] = animationJson;

	std::ostringstream stringStream;
	stringStream << newName << ".bonemf";

	//std::ofstream o(fileName);

	if (!doJson) {
		std::string fileName = stringStream.str();
		std::ofstream o(fileName, std::ios::out | std::ios::binary);
		std::vector<std::uint8_t> v_cbor = json::to_cbor(modelFile);
		o.write((char*)&v_cbor[0], v_cbor.size() * sizeof(std::uint8_t));
		o.close();
	}
	else {
		stringStream << ".json";
		std::string fileName = stringStream.str();
		std::ofstream o(fileName);
		o << std::setw(4) << modelFile << std::endl;
		o.close();
	}


	lSdkManager->Destroy();
	return 0;
}
