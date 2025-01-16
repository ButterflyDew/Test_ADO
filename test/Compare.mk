CXX = g++ # C++编译器
CXXFLAGS = -std=c++20 -O2 -g# 编译选项

TARGET = CompareTest # 目标可执行文件名

# 列出你的源文件
SOURCES = ../AdjacencyList.hpp ../HL/HL_basic.hpp ../ADO/ADO_basic.hpp ../ADO/ADO_fixed_Ak.hpp ../Utilities.hpp Load_Graph.hpp Compare.cpp

# 生成目标的规则
$(TARGET): $(SOURCES) $(HEADERS)
	$(CXX) $(CXXFLAGS) -o $@ $(SOURCES)

.PHONY: clean

clean:
	rm -f $(TARGET)
