CXX = g++ # C++编译器
CXXFLAGS = -std=c++17 -O2# 编译选项

TARGET = HLTest # 目标可执行文件名

# 列出你的源文件
SOURCES = ../AdjacencyList.hpp ../HL/HL_basic.hpp ../Utilities.hpp Load_Graph.hpp HLReadTest.cpp

# 生成目标的规则
$(TARGET): $(SOURCES) $(HEADERS)
	$(CXX) $(CXXFLAGS) -o $@ $(SOURCES)

.PHONY: clean

clean:
	rm -f $(TARGET)
