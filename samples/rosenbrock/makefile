# Copyright (C) 2013 Gilberto Noronha
#
# This is free software; you can redistribute it and/or modify it under the
# terms of the GNU General Public License as published by the Free Software
# Foundation; either version 3 of the License, or (at your option) any later
# version.

# This software is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
# General Public License for more details.

# You should have received a copy of the GNU General Public License along
# with this program; if not, write to the Free Software Foundation, Inc., 51
# Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.

INCLUDE_PATH = ../../src
CFLAGS = -I$(INCLUDE_PATH) -Wall -Wextra -pedantic -O2 -std=c99

TARGET = main.out

all: $(TARGET)

$(TARGET): main.o
	$(CC) $(CFLAGS) main.o -o main.out

main.o: main.c
	$(CC) $(CFLAGS) -c main.c -o main.o 

.PHONY: all clean

clean:
	rm -f main.o main.out
