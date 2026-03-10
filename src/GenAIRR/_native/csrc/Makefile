SHELL := /bin/bash

BUILD_DIR := build

.PHONY: all test clean rebuild

all:
	@mkdir -p $(BUILD_DIR)
	@cd $(BUILD_DIR) && cmake .. -DCMAKE_BUILD_TYPE=Debug -Wno-dev 2>/dev/null
	@$(MAKE) -C $(BUILD_DIR) --no-print-directory

test: all
	@cd $(BUILD_DIR) && ctest --output-on-failure

clean:
	@rm -rf $(BUILD_DIR)

rebuild: clean all
