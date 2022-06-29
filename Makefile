TARGET = smallScaleTest

all:
	$(MAKE) clean
	$(MAKE) -C source
clean:
	$(RM) $(TARGET) ./source/$(TARGET)
